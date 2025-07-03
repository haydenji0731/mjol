from mjol.base import *
import pandas as pd
from typing import List

HDR = [
    'chr', 'src', 'feature_type', 'start', 
    'end', 'score', 'strand', 'frame', 'attributes'
]
FEATURES = ['gene', 'transcript', 'exon', 'cds']

class GAn(BaseModel):
    filename : str
    format : str
    genes : Optional[Dict[str, Dict[str, GFeature]]] = None
    txes : Optional[Dict[str, Dict[str, GFeature]]] = None
    orphans : Optional[List[GFeature]] = None

    def build_db(self):
        in_df = pd.read_csv(self.filename, sep='\t', 
                        comment='#', header=None)
        in_df.columns = HDR

        feature_types = in_df['feature_type'].unique()

        # TODO: better logic?
        # if feature_types != FEATURES:
        #     raise ValueError(f"unexpected feature types: {feature_types}")

        sub_dfs = dict()
        for x in FEATURES:
            sub_dfs[x] = in_df[in_df['feature_type'] == x]

        gene_df = sub_dfs['gene']
        for chr, curr_df in gene_df.groupby("chr"):
            self.genes[chr] = self._pd_df2genes(curr_df)
        
        tx_df = sub_dfs['transcript']
        for chr, curr_df in tx_df.groupby("chr"):
            self.txes[chr] = self._pd_df2txes(curr_df)
    
    def _pd_df2genes(self, df):
        genes = dict()
        for _, row in df.iterrows():
            gfeat = self._pd_row2gfeature(row)
            genes[gfeat.id] = gfeat
        return genes
    
    def _pd_df2txes(self, chr, df):
        txes = dict()
        for _, row in df.iterrows():
            gfeat = self._pd_row2gfeature(row)
            if not gfeat.parent:
                self.orphans.append(gfeat)
            self.genes[chr][gfeat.parent].children[gfeat.id] = gfeat
            txes[gfeat.id] = gfeat
        return txes
    
    def _pd_df2exons(self, chr, df):
        raise NotImplementedError
    
    def _pd_df2cdses(self, chr, df):
        raise NotImplementedError
    
    def _pd_row2gfeature(self, row) -> GFeature:
        try:
            gfeat = GFeature(
                chr = row['chr'],
                src = row['src'],
                feature_type = row['feature_type'],
                start = row['start'],
                end = row['end'],
                score = None if row['score'] == '.' else row['score'],
                strand = row['strand'],
                frame = row['frame'],
                attributes = load_attributes(row['attributes'], 
                                        kv_sep = ' ' if self.format == '' else '=')
            )
        except Exception as e:
            raise RuntimeError(f"error while converting pd series to GFeature: {e}")
        return gfeat