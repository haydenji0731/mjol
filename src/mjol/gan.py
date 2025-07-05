from mjol.base import *
import pandas as pd
from typing import List
from concurrent.futures import ProcessPoolExecutor

HDR = [
    'chr', 'src', 'feature_type', 'start', 
    'end', 'score', 'strand', 'frame', 'attributes'
]

# TODO: make this more flexible
FEATURES = ['gene', 'transcript', 'exon', 'CDS']

class GAn(BaseModel):
    filename : str
    format : str
    genes : Dict[str, Dict[str, GFeature]] = Field(default_factory=dict)
    txes : Dict[str, Dict[str, GFeature]] = Field(default_factory=dict)
    orphans : List[GFeature] = Field(default_factory=list)

    def build_db(self, n_threads : int = 1):
        in_df = pd.read_csv(self.filename, sep='\t', 
                        comment='#', header=None)
        in_df.columns = HDR

        sub_dfs = dict()
        for x in FEATURES:
            sub_dfs[x] = in_df[in_df['feature_type'] == x]

        # process genes
        gene_df = sub_dfs['gene']
        for chr, curr_df in gene_df.groupby("chr"):
            self.genes[chr] = self._pd_df2genes(curr_df)
        
        # process transcripts
        tx_df = sub_dfs['transcript']
        for chr, curr_df in tx_df.groupby("chr"):
            self.txes[chr] = self._pd_df2txes(curr_df)
        
        # process exons
        exon_df = sub_dfs['exon']
        dfs_by_chr = [x for _, x in exon_df.groupby("chr")]
        chrs = [x for x, _ in exon_df.groupby("chr")]
        if n_threads > 1:
            with ProcessPoolExecutor(max_workers=n_threads) as executor:
                results = list(executor.map(self._pd_df2exons_pll, dfs_by_chr))
            for i in range(len(chrs)):
                orphans, exons = results[i]
                self.orphans += orphans
                for id in exons:
                    self.txes[chrs[i]][id].children += exons[id]

        else:
            for chr, curr_df in exon_df.groupby("chr"):
                self._pd_df2exons(curr_df)

        # process cdses
        cds_df = sub_dfs['CDS']
        
                
    
    def _pd_df2genes(self, df):
        genes = dict()
        for _, row in df.iterrows():
            gfeat = self._pd_row2gfeature(row)
            genes[gfeat.id] = gfeat
        return genes
    
    def _pd_df2txes(self, df):
        txes = dict()
        for _, row in df.iterrows():
            gfeat = self._pd_row2gfeature(row)
            if not gfeat.parent:
                self.orphans.append(gfeat)
                continue
            self.genes[gfeat.chr][gfeat.parent].children.append(gfeat)
            txes[gfeat.id] = gfeat
        return txes
    
    def _pd_df2exons(self, df):
        for _, row in df.iterrows():
            gfeat = self._pd_row2gfeature(row)
            if not gfeat.parent:
                self.orphans.append(gfeat)
                continue
            self.txes[gfeat.chr][gfeat.parent].children.append(gfeat)
    
    def _pd_df2exons_pll(self, df):
        orphans = []
        exons = dict()
        for _, row in df.iterrows():
            gfeat = self._pd_row2gfeature(row)
            if not gfeat.parent:
                orphans.append(gfeat)
                continue
            if gfeat.parent not in exons:
                exons[gfeat.parent] = [gfeat]
            else:
                exons[gfeat.parent].append(gfeat)
        return orphans, exons

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
                                        kv_sep = ' ' if self.format == '' else '='),
                children = []
            )
        except Exception as e:
            raise RuntimeError(f"error while converting pd series to GFeature: {e}")
        return gfeat
    
    # TODO: add option to sort
    def to_gff(self, fp: str):
        with open(fp, "w") as f:
            f.write("##gff-version 3\n")
            for gene_chr in self.genes.values():
                for gene in gene_chr.values():
                    f.write(gene.to_gff_entry())
                    for transcript in gene.children or []:
                        f.write(transcript.to_gff_entry())
                        for exon in transcript.children or []:
                            f.write(exon.to_gff_entry())
                            # CDS not implemented yet
                            for cds in exon.children or []:
                                f.write(cds.to_gff_entry())
            for orphan in self.orphans:
                f.write(orphan.to_gff_entry())

