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

class GAn(BaseModel):
    filename : str
    format : str



    
    features : List[str] = ['gene', 'transcript', 'exon', 'CDS']
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
        for x in self.features:
            sub_dfs[x] = in_df[in_df['feature_type'] == x]
            
        gene_df = sub_dfs[self.features[0]]
        tx_df = sub_dfs[self.features[1]]
        exon_df = sub_dfs[self.features[2]]
        cds_df = sub_dfs[self.features[3]]

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

    def delete_gene(self, id:str):
        for _, genes_dict in self.genes.items():
            if id in genes_dict:
                delete_gene = genes_dict[id]
                for child in delete_gene.children:
                    del self.txes[delete_gene.chr][child.id]
                del genes_dict[id]
                return delete_gene.to_gff_entry(children=True)
        raise KeyError(f'Gene ID {id} not found in the gene annotation')

    def add_gene(self, gene:GFeature) -> str:
        for chr in self.genes.keys():
            if gene.id in self.genes[chr]:
                print("WARNING: A gene with the same ID exists and will be overwritten")
        self.genes[gene.chr][gene.id] = gene
        for child in gene.children:
            self.txes[child.chr][child.id] = child
        return gene.to_gff_entry(children= True)


    def merge(self, other:GAnRef):

        # check for duplicate ids
        self_genes = set()
        self_txes = set()
        other_genes = set()
        other_txes = set()
        for chromosome in self.genes.keys():
            self_genes.union(set(self.genes[chromosome].keys()))
            self_txes.union(set(self.txes[chromosome].keys()))
        for chromosome in other.genes.keys():
            other_genes.union(set(other.genes[chromosome].keys()))
            other_txes.union(set(other.txes[chromosome].keys()))
        if (self_genes & other_genes) or (self_txes & other_txes):
            print(f"WARNING: {len(self_genes & other_genes) + len(self_txes & other_txes)} duplicate IDs found")

        for chromosome in (self.genes.keys() | other.genes.keys()):
            self.genes[chromosome] = {**self.genes[chromosome], **other.genes[chromosome]}
        for chromosome in (self.txes.keys() | other.txes.keys()):
            self.txes[chromosome] = {**self.txes[chromosome], **other.txes[chromosome]}
        self.filename = self.filename + other.filename
        return self
    
    # returns a reference
    def get_gene(self, id:str):
        for chr, genes_dict in self.genes.items():
            if id in genes_dict:
                return genes_dict[id]
        raise KeyError(f'Gene ID {id} not found in the gene annotation')
    
    def _get_chr_gene(self, id:str) -> tuple:
        for chr, genes_dict in self.genes.items():
            if id in genes_dict:
                return chr, genes_dict[id]
        raise KeyError(f'Gene ID {id} not found in the gene annotation')
    
    def get_descendants(self, feature: GFeature) -> List[GFeature]:
        descendants = list()
        for child in feature.children:
            descendants.append(child)
            descendants.extend(self.get_descendants(child))
        return descendants
        
    def save_as_gix(self, filepath : str):
        with open(filepath, 'wb') as fh:
            pickle.dump(self, fh)
    
    # NOTE: this function can only change IDs for genes
    def solve_synonym(self, self_id: str, other: GAnRef, other_id: str,
                      update_gene_attributes: List[tuple] = [], update_children_attributes: List[tuple] = [],
                      exclude_attributes: List[str] = []) -> tuple:
        # update gene GFeature
        old_chr, old_feature = self._get_chr_gene(self_id)
        new_chr, new_feature = other._get_chr_gene(other_id)

        # save old entries
        old_entries = old_feature.to_gff_entry(children=True)
        # gene-level updates
        ## class attribute

        old_feature.id = new_feature.id
        for key_old, key_new in update_gene_attributes:
            if (key_old in old_feature.attributes) and (key_new in new_feature.attributes):
                old_feature.attributes[key_old] = new_feature.attributes[key_new]
        for key in exclude_attributes:
            if key in old_feature.attributes:
                del old_feature.attributes[key]
            
        ## gene dictionary 

        del self.genes[old_chr][self_id]
        self.genes[old_chr][other_id] = old_feature
        ## change parent field of transcripts

        for child in old_feature.children:
            ## class attribute
            child.parent = new_feature.id
            child.attributes['Parent'] = new_feature.id
        # update attributes for children
        for descendant in self.get_descendants(old_feature):
            for key_old, key_new in update_children_attributes:
                if (key_old in descendant.attributes) and (key_new in new_feature.attributes):
                    descendant.attributes[key_old] = new_feature.attributes[key_new]
            for key in exclude_attributes:
                if key in descendant.attributes:
                    del descendant.attributes[key]
        # get updated entries
        new_entries = old_feature.to_gff_entry(children=True)

        return (old_entries, new_entries)
        
def load_gan_from_gix(filepath : str):
    try:
        with open(filepath, 'rb') as fh:
            res = pickle.load(fh)
    except Exception as e:
        raise RuntimeError(f"error while loading .gix: {e}")
    return res


