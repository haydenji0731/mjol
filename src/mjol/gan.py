from mjol.base import *
import pandas as pd
from typing import List
from concurrent.futures import ProcessPoolExecutor
import numpy as np

GAnRef = ForwardRef("GAn")

HDR = [
    'chr', 'src', 'feature_type', 'start', 
    'end', 'score', 'strand', 'frame', 'attributes'
]

# TODO: make this more flexible
FEATURES = ['gene', 'transcript', 'exon', 'CDS']

# make entries in a column unique (for IDs)
def make_unique(df: pd.DataFrame , column: str) -> tuple[pd.DataFrame, np.ndarray]:
    duplicates = df.loc[df[column].duplicated(), column].unique()

    count_series = df.groupby(column).cumcount()
    new_id = df[column] + np.where(count_series > 0, "-" + count_series.astype(str), '')
    df.loc[:, column] = new_id

    return (df, duplicates)

class GAn(BaseModel):
    filename : str
    format : str
    genes : Dict[str, Dict[str, GFeature]] = Field(default_factory=dict)
    txes : Dict[str, Dict[str, GFeature]] = Field(default_factory=dict)
    orphans : List[GFeature] = Field(default_factory=list)

    # duplicate_id: 'make_unique' or 'raise_error'
    def build_db(self, duplicate_id : str = 'make_unique', n_threads : int = 1):
        in_df = pd.read_csv(self.filename, sep='\t', 
                        comment='#', header=None)
        in_df.columns = HDR
        
        # parse attributes
        in_df['attributes'] = in_df['attributes'].apply(
            lambda s: load_attributes(s, kv_sep = ' ' if self.format == '' else '=')
        )
        in_df['ID'] = in_df['attributes'].apply(lambda s: s["ID"] if "ID" in s.keys() else None)

        sub_dfs = dict()
        for x in FEATURES:
            sub_dfs[x] = in_df[in_df['feature_type'] == x]
            
        # process genes and transcripts
        gene_df = sub_dfs['gene']
        tx_df = sub_dfs['transcript']
    
        # check ids
        if duplicate_id == 'raise_error':
            _, gene_duplicates = make_unique(gene_df, "ID")
            _, tx_duplicates = make_unique(tx_df, "ID")
            if gene_duplicates or tx_duplicates:
                raise ValueError("There are duplicate IDs. Set `duplicate_id` to 'make_unique' to avoid this" )
        elif duplicate_id == 'make_unique':
            gene_df, _ = make_unique(gene_df, "ID")
            tx_df, _ = make_unique(tx_df, "ID")
        else:
            raise ValueError(f"duplicate_id must be 'raise_error' or 'make_unique', got {duplicate_id}")
        
        gene_df.loc[:, 'attributes'] = gene_df.apply(
            lambda row: {**row['attributes'], "ID": row["ID"]},
            axis=1
        )
        tx_df.loc[:,'attributes'] = tx_df.apply(
            lambda row: {**row['attributes'], 'ID':row['ID']},
            axis=1
        )
        # add genes to dict
        for chr, curr_df in gene_df.groupby("chr"):
            self.genes[chr] = self._pd_df2genes(curr_df)
        # add transcripts to dict
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
                attributes = row['attributes'], 
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
                    f.write(gene.to_gff_entry(children=True))
            for orphan in self.orphans:
                f.write(orphan.to_gff_entry())

    # TODO: transcripts might still be stored; it won't show up via to_gff() but still exists in memory
    def delete_gene(self, id:str):
        for _, genes_dict in self.genes.items():
            if id in genes_dict:
                delete_gene = genes_dict[id]
                for child in delete_gene.children:
                    del self.txes[child.id]
                del genes_dict[id]
                return delete_gene.to_gff_entry(children=True)
        raise KeyError(f'Gene ID {id} not found in the gene annotation')

    # TODO: make_unique option to cover cases where ids may not be unique
    def merge(self, other: GAnRef, duplicate_id: str = 'make_unique', make_unique_suffix: str = ('_1')) -> GAnRef:
        for chromosome in (self.genes.keys() | other.genes.keys()):
            if (chromosome in self.genes.keys()) and (chromosome in other.genes.keys()):
                overlapping_ids = self.genes[chromosome].keys() & other.genes[chromosome].keys()
                if overlapping_ids:
                    if duplicate_id == 'raise_error':
                        raise ValueError(f"Duplicate gene IDs found on chromosome {chromosome}: {overlapping_ids}")
                    elif duplicate_id == 'make_unique':
                        temp_dict = {
                            (key + make_unique_suffix if key in overlapping_ids else key): value
                            for key, value in other.genes[chromosome].items()
                        }
                        self.genes[chromosome] = {**self.genes[chromosome], **temp_dict}
                    else:
                        raise ValueError(f"duplicate_id must be 'raise_error' or 'make_unique', got {duplicate_id}")
                else:
                    self.genes[chromosome] = {**self.genes[chromosome], **other.genes[chromosome]}
        for chromosome in (self.txes.keys() | other.txes.keys()):
            if (chromosome in self.txes.keys()) and (chromosome in other.txes.keys()):
                overlapping_ids = self.txes[chromosome].keys() & other.txes[chromosome].keys()
                if overlapping_ids:
                    if duplicate_id == 'raise_error':
                        raise ValueError(f"Duplicate gene IDs found on chromosome {chromosome}: {overlapping_ids}")
                    elif duplicate_id == 'make_unique':
                        temp_dict = {
                            (key + make_unique_suffix if key in overlapping_ids else key): value
                            for key, value in other.txes[chromosome].items()
                        }
                        self.txes[chromosome] = {**self.txes[chromosome], **temp_dict}
                    else:
                        raise ValueError(f"duplicate_id must be 'raise_error' or 'make_unique', got {duplicate_id}")
                else:
                    self.txes[chromosome] = {**self.txes[chromosome], **other.txes[chromosome]}
        return self

