from mjol.base import *
from mjol.utils import *
import pandas as pd

HDR = [
    'chr', 'src', 'feature_type', 'start', 
    'end', 'score', 'strand', 'frame', 'attributes'
]

class GAn(BaseModel):
    file_name : str
    file_fmt : str
    iak : str = 'id'
    pak : str = 'parent'
    features : dict = Field(default_factory=dict)
    lookup : dict = Field(default_factory=dict)
    
    def build_db(self):
        in_df = pd.read_csv(self.file_name, sep='\t', comment='#', header=None)
        in_df.columns = HDR
        in_df['attributes'] = in_df['attributes'].apply(
            lambda s : load_attributes(
                s, kv_sep=' ' if self.file_fmt.lower() == 'gtf' else '='
            )
        )
        
        rows = in_df.to_dict('records')
        for row in rows:
            f = self._create_gfeature(row)
            self.features[f.uid] = f
            if f.aid:

                if f.aid in self.lookup:
                    self.lookup[f.aid].append(f.uid)

                self.lookup[f.aid] = [f.uid]

            if f.parent in self.lookup:
                puid = self.get_uid(f.parent, f)
                
                f.set_parent_uid(puid)
                self.features[puid].add_a_child(f)
    
    # collision-safe
    def get_uid(self, aid : str, f : GFeature = None):
        uids = self.lookup[aid]
        if len(uids) == 1:
            return uids[0]
        if not f:
            raise RuntimeError(f'provide a feature to resolve lookup collision')
        return max(uids, key=lambda uid: f.calc_sim(self.features[uid]))
        
    def get_feature(self, uid : str):
        if uid not in self.features:
            raise KeyError(f'{uid} not found in features')
        return self.features[uid]
    
    def get_desc(self, uid : str):
        res = []
        for child in self.features[uid].children:
            res.append(child)
            res.extend(self.get_desc(child))
        return res
    
    # TODO
    def delete_feature(self):
        raise NotImplementedError

    # TODO
    def add_feature(self):
        raise NotImplementedError
    
    def save_as_gix(self, file_path : str):
        with open(file_path, 'wb') as fh:
            pickle.dump(self, fh)

    def _create_gfeature(self, row):
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
                    iak = self.iak,
                    pak = self.pak
                )
        return gfeat

def load_from_gix(file_path : str):
    try:
        with open(file_path, 'rb') as fh:
            res = pickle.load(fh)
    except Exception as e:
        raise RuntimeError(f"error while loading {file_path} : {e}")
    return res