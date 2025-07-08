from mjol.base import *
import pandas as pd
from typing import List
from concurrent.futures import ProcessPoolExecutor
import pickle

HDR = [
    'chr', 'src', 'feature_type', 'start', 
    'end', 'score', 'strand', 'frame', 'attributes'
]

# TODO: abstractify
FEATURES = ['gene', 'transcript', 'exon', 'CDS']

class GAn(BaseModel):
    filename : str
    format : str
    genes : Dict[str, Dict[str, GFeature]] = Field(default_factory=dict)
    txes : Dict[str, Dict[str, GFeature]] = Field(default_factory=dict)
    orphans : List[GFeature] = Field(default_factory=list)

    def build_db(self, use_iv : bool = False, n_threads : int = 1):
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
            self.txes[chr] = self._pd_df2txes(curr_df, use_iv = use_iv)
        
        print("exon")
        # process exons
        exon_df = sub_dfs['exon']

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

        print("cds")
        cds_df = sub_dfs['CDS']

        if use_iv:
            # TODO: why is this slow?
            # if n_threads > 1:
            #     chrs, dfs_by_chr = zip(*exon_df.groupby("chr"))
            #     with ProcessPoolExecutor(max_workers=n_threads) as executor:
            #         results = list(executor.map(self._pd_df2cdses_pll, dfs_by_chr))
            #     for i in range(len(chrs)):
            #         orphans, cdses = results[i]
            #         self.orphans += orphans
            #         for tx_id in cdses:
            #             for cx in cdses[tx_id]:
            #                 parent_exon = self.txes[cx.chr][cx.parent].query_itree(cx.start, cx.end)
            #                 try:
            #                     assert len(parent_exon) == 0 # sanity check
            #                 except Exception as e:
            #                     print(parent_exon)
            #                     raise ValueError
            #                 next(iter(parent_exon)).add_a_child(cx)
            # else:
            self._pd_df2cdses(cds_df, use_iv = True)
        else:
            # TODO: why is this slow?
            # if n_threads > 1:
            #     chrs, dfs_by_chr = zip(*exon_df.groupby("chr"))
            #     with ProcessPoolExecutor(max_workers=n_threads) as executor:
            #         results = list(executor.map(self._pd_df2cdses_pll, dfs_by_chr))
            #     for i in range(len(chrs)):
            #         orphans, cdses = results[i]
            #         self.orphans += orphans
            #         for tx_id in cdses:
            #             for cx in cdses[tx_id]:
            #                 self.txes[cx.chr][cx.parent].add_a_child(cx)
            # else:
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

    def _pd_df2cdses_pll(self, df):
        orphans = []
        cdses = dict()
        for _, row in df.iterrows():
            cds = self._pd_row2gfeature(row)
            if not cds.parent:
                orphans.append(self)
                continue
            if cds.parent not in cdses:
                cdses[cds.parent] = [cds]
            else:
                cdses[cds.parent].append(cds)
        return orphans, cdses
    
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
                    attributes = load_attributes(row['attributes'], 
                                            kv_sep = ' ' if self.format == '' else '='),
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
    
    def save_as_gix(self, filepath : str):
        with open(filepath, 'wb') as fh:
            pickle.dump(self, fh)
    
def load_gan_from_gix(filepath : str):
    try:
        with open(filepath, 'rb') as fh:
            res = pickle.load(fh)
    except Exception as e:
        raise RuntimeError(f"error while loading .gix: {e}")
