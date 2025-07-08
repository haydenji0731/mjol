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