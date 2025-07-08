from mjol.base import *
import pandas as pd
from typing import List
from concurrent.futures import ProcessPoolExecutor
import pickle
import numpy as np

GAnRef = ForwardRef("GAn")

HDR = [
    'chr', 'src', 'feature_type', 'start', 
    'end', 'score', 'strand', 'frame', 'attributes'
]

FEATURES = ['gene', 'transcript', 'exon', 'CDS']

class GAn(BaseModel):
    filename : str
    format : str

    genes : Dict[str, Dict[str, GFeature]] = Field(default_factory=dict)
    txes : Dict[str, Dict[str, GFeature]] = Field(default_factory=dict)
    orphans : List[GFeature] = Field(default_factory=list)

    # duplicate_id: 'make_unique' or 'raise_error'
    def build_db(self, use_iv : bool = False, n_threads : int = 1):
        in_df = pd.read_csv(self.filename, sep='\t', 
                        comment='#', header=None)
        in_df.columns = HDR
        
        # parse attributes
        in_df['attributes'] = in_df['attributes'].apply(
            lambda s: load_attributes(s, kv_sep = ' ' if self.format == '' else '=')
        )

        sub_dfs = dict()
        for x in FEATURES:
            sub_dfs[x] = in_df[in_df['feature_type'] == x]
            
        gene_df = sub_dfs['gene']
        tx_df = sub_dfs['transcript']
        exon_df = sub_dfs['exon']
        cds_df = sub_dfs['CDS']

        # process genes
        for chr, curr_df in gene_df.groupby("chr"):
            self.genes[chr] = self._pd_df2genes(curr_df)

        # process transcripts
        for chr, curr_df in tx_df.groupby("chr"):
            self.txes[chr] = self._pd_df2txes(curr_df, use_iv = use_iv)
        
        # process exons
        if n_threads > 1:
            chrs, dfs_by_chr = zip(*exon_df.groupby("chr"))
            with ProcessPoolExecutor(max_workers=n_threads) as executor:
                results = list(executor.map(self._pd_df2exons_pll, dfs_by_chr))
            for i in range(len(chrs)):
                orphans, exons = results[i]
                self.orphans += orphans
                for tx_id in exons:
                    for ex in exons[tx_id]:
                        self.txes[chrs[i]][tx_id].add_a_child(ex)

        else:
            for chr, curr_df in exon_df.groupby("chr"):
                self._pd_df2exons(curr_df)

        # process cdses
        if use_iv:
            self._pd_df2cdses(cds_df, use_iv = True)
        else:
            self._pd_df2cdses(cds_df)
    
    def _pd_df2genes(self, df):
        genes = dict()
        for _, row in df.iterrows():
            gene = self._pd_row2gfeature(row)
            genes[gene.id] = gene
        return genes
    
    def _pd_df2txes(self, df, use_iv : bool = False):
        txes = dict()
        for _, row in df.iterrows():
            tx = self._pd_row2gfeature(row, use_iv = use_iv)
            if not tx.parent:
                self.orphans.append(tx)
                continue
            self.genes[tx.chr][tx.parent].add_a_child(tx)
            txes[tx.id] = tx
        return txes
    
    def _pd_df2exons(self, df):
        for _, row in df.iterrows():
            exon = self._pd_row2gfeature(row)
            if not exon.parent:
                self.orphans.append(exon)
                continue
            self.txes[exon.chr][exon.parent].add_a_child(exon)
    
    def _pd_df2exons_pll(self, df):
        orphans = []
        exons = dict()
        for _, row in df.iterrows():
            exon = self._pd_row2gfeature(row)
            if not exon.parent:
                orphans.append(self)
                continue
            if exon.parent not in exons:
                exons[exon.parent] = [exon]
            else:
                exons[exon.parent].append(exon)
        return orphans, exons
    
    def _pd_df2cdses(self, df, use_iv : bool = False):
        for _, row in df.iterrows():
            cds = self._pd_row2gfeature(row)
            if not cds.parent:
                self.orphans.append(cds)
                continue
            if use_iv:
                parent_exon = self.txes[cds.chr][cds.parent].query_itree(cds.start, cds.end)
                try:
                    assert len(parent_exon) == 0 # sanity check
                except Exception as e:
                    print(parent_exon)
                    raise ValueError
                
                next(iter(parent_exon)).add_a_child(cds)
            else:
                # TODO: update divider
                self.txes[cds.chr][cds.parent].add_a_child(cds)
    
    def _pd_row2gfeature(self, row, use_iv : bool = False) -> GFeature:
        try:
            if use_iv:
                gfeat = IntervalGFeature(
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
            else:
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
                    del self.txes[delete_gene.chr][child.id]
                del genes_dict[id]
                return delete_gene.to_gff_entry(children=True)
        raise KeyError(f'Gene ID {id} not found in the gene annotation')

    # NOTE: if self and other have entities with same ID but different loc, problem
    def merge(self, other):
        # TODO: compute the intersection between self and other's key space
        for chromosome in (self.genes.keys() | other.genes.keys()):
            self.genes[chromosome] = {**self.genes[chromosome], **other.genes[chromosome]}
        for chromosome in (self.txes.keys() | other.txes.keys()):
            self.txes[chromosome] = {**self.txes[chromosome], **other.txes[chromosome]}
        return self
    
    def save_as_gix(self, filepath : str):
        with open(filepath, 'wb') as fh:
            pickle.dump(self, fh)
    
def load_gan_from_gix(filepath : str):
    try:
        with open(filepath, 'rb') as fh:
            res = pickle.load(fh)
    except Exception as e:
        raise RuntimeError(f"error while loading .gix: {e}")
    return res