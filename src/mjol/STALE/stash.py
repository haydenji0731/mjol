# make entries in a column unique (for IDs)
# def make_unique(df: pd.DataFrame , column: str) -> tuple[pd.DataFrame, np.ndarray]:
#     duplicates = df.loc[df[column].duplicated(), column].unique()

#     count_series = df.groupby(column).cumcount()
#     new_id = df[column] + np.where(count_series > 0, "-" + count_series.astype(str), '')
#     df.loc[:, column] = new_id

#     return (df, duplicates)

# def _pd_df2cdses_pll(self, df):
#     orphans = []
#     cdses = dict()
#     for _, row in df.iterrows():
#         cds = self._pd_row2gfeature(row)
#         if not cds.parent:
#             orphans.append(self)
#             continue
#         if cds.parent not in cdses:
#             cdses[cds.parent] = [cds]
#         else:
#             cdses[cds.parent].append(cds)
#     return orphans, cdses

# class GAn(BaseModel):
#     filename : str
#     format : str
#     genes : Dict[str, Dict[str, GFeature]] = Field(default_factory=dict)
#     txes : Dict[str, Dict[str, GFeature]] = Field(default_factory=dict)
#     orphans : List[GFeature] = Field(default_factory=list)

#     # duplicate_id: 'make_unique' or 'raise_error'
#     def build_db(self, use_iv : bool = False, n_threads : int = 1):
#         in_df = pd.read_csv(self.filename, sep='\t', 
#                         comment='#', header=None)
#         in_df.columns = HDR
        
#         # parse attributes
#         in_df['attributes'] = in_df['attributes'].apply(
#             lambda s: load_attributes(s, kv_sep = ' ' if self.format == '' else '=')
#         )
#         # in_df['ID'] = in_df['attributes'].apply(lambda s: s["ID"] if "ID" in s.keys() else None)

#         sub_dfs = dict()
#         for x in FEATURES:
#             sub_dfs[x] = in_df[in_df['feature_type'] == x]
            
#         # process genes and transcripts
#         gene_df = sub_dfs['gene']
#         tx_df = sub_dfs['transcript']
    
#         # check ids
#         # if duplicate_id == 'raise_error':
#         #     _, gene_duplicates = make_unique(gene_df, "ID")
#         #     _, tx_duplicates = make_unique(tx_df, "ID")
#         #     if gene_duplicates or tx_duplicates:
#         #         raise ValueError("There are duplicate IDs. Set `duplicate_id` to 'make_unique' to avoid this" )
#         # elif duplicate_id == 'make_unique':
#         #     gene_df, _ = make_unique(gene_df, "ID")
#         #     tx_df, _ = make_unique(tx_df, "ID")
#         # else:
#         #     raise ValueError(f"duplicate_id must be 'raise_error' or 'make_unique', got {duplicate_id}")
        
#         # gene_df.loc[:, 'attributes'] = gene_df.apply(
#         #     lambda row: {**row['attributes'], "ID": row["ID"]},
#         #     axis=1
#         # )
#         # tx_df.loc[:,'attributes'] = tx_df.apply(
#         #     lambda row: {**row['attributes'], 'ID':row['ID']},
#         #     axis=1
#         # )
#         # add genes to dict
#         for chr, curr_df in gene_df.groupby("chr"):
#             self.genes[chr] = self._pd_df2genes(curr_df)
#         # add transcripts to dict
#         tx_df = sub_dfs['transcript']
#         for chr, curr_df in tx_df.groupby("chr"):
#             self.txes[chr] = self._pd_df2txes(curr_df, use_iv = use_iv)
        
#         print("exon")
#         # process exons
#         exon_df = sub_dfs['exon']

#         if n_threads > 1:
#             chrs, dfs_by_chr = zip(*exon_df.groupby("chr"))
#             with ProcessPoolExecutor(max_workers=n_threads) as executor:
#                 results = list(executor.map(self._pd_df2exons_pll, dfs_by_chr))
#             for i in range(len(chrs)):
#                 orphans, exons = results[i]
#                 self.orphans += orphans
#                 for tx_id in exons:
#                     for ex in exons[tx_id]:
#                         self.txes[chrs[i]][tx_id].add_a_child(ex)

#         else:
#             for chr, curr_df in exon_df.groupby("chr"):
#                 self._pd_df2exons(curr_df)

#         print("cds")
#         cds_df = sub_dfs['CDS']

#         if use_iv:
#             # TODO: why is this slow?
#             # if n_threads > 1:
#             #     chrs, dfs_by_chr = zip(*exon_df.groupby("chr"))
#             #     with ProcessPoolExecutor(max_workers=n_threads) as executor:
#             #         results = list(executor.map(self._pd_df2cdses_pll, dfs_by_chr))
#             #     for i in range(len(chrs)):
#             #         orphans, cdses = results[i]
#             #         self.orphans += orphans
#             #         for tx_id in cdses:
#             #             for cx in cdses[tx_id]:
#             #                 parent_exon = self.txes[cx.chr][cx.parent].query_itree(cx.start, cx.end)
#             #                 try:
#             #                     assert len(parent_exon) == 0 # sanity check
#             #                 except Exception as e:
#             #                     print(parent_exon)
#             #                     raise ValueError
#             #                 next(iter(parent_exon)).add_a_child(cx)
#             # else:
#             self._pd_df2cdses(cds_df, use_iv = True)
#         else:
#             # TODO: why is this slow?
#             # if n_threads > 1:
#             #     chrs, dfs_by_chr = zip(*exon_df.groupby("chr"))
#             #     with ProcessPoolExecutor(max_workers=n_threads) as executor:
#             #         results = list(executor.map(self._pd_df2cdses_pll, dfs_by_chr))
#             #     for i in range(len(chrs)):
#             #         orphans, cdses = results[i]
#             #         self.orphans += orphans
#             #         for tx_id in cdses:
#             #             for cx in cdses[tx_id]:
#             #                 self.txes[cx.chr][cx.parent].add_a_child(cx)
#             # else:
#             self._pd_df2cdses(cds_df)

            # if (chromosome in self.txes.keys()) and (chromosome in other.txes.keys()):
            #     overlapping_ids = self.txes[chromosome].keys() & other.txes[chromosome].keys()
            #     if overlapping_ids:
            #         if duplicate_id == 'raise_error':
            #             raise ValueError(f"Duplicate gene IDs found on chromosome {chromosome}: {overlapping_ids}")
            #         elif duplicate_id == 'make_unique':
            #             temp_dict = {
            #                 (key + make_unique_suffix if key in overlapping_ids else key): value
            #                 for key, value in other.txes[chromosome].items()
            #             }
            #             self.txes[chromosome] = {**self.txes[chromosome], **temp_dict}
            #         else:
            #             raise ValueError(f"duplicate_id must be 'raise_error' or 'make_unique', got {duplicate_id}")
            #     else: