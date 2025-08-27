from mjol.base import *
from mjol.utils import *
import pandas as pd
import pickle

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
    
    # NOTE: child feature must come after parent feature in GFF file
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
            if f.uid in self.features:
                raise RuntimeError(f'non-unique uid detected : {f.uid}')
            
            self.features[f.uid] = f

            if f.aid:
                # attribute ID collision detected
                if f.aid in self.lookup:
                    self.lookup[f.aid].append(f.uid)
                else:
                    self.lookup[f.aid] = [f.uid]
        
        for f in self.features.values():
            if f.paid in self.lookup:
                puid = self.get_uid(f.paid, f)
                f.set_parent_uid(puid)
                parent = self.get_feature(puid)
                parent.add_a_child(f)

            # OLD version for constructing parent-child relationship
            # if f.parent in self.lookup:
            #     puid = self.get_uid(f.parent, f)
            #     f.set_parent_uid(puid)
            #     parent = self.get_feature(puid)
            #     parent.add_a_child(f)
    
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
    
    def pop_feature(self, uid: str, include_children = True) -> str:
        if uid not in self.features:
            raise KeyError(f'{uid} not found in features')
        delete_feature = self.features[uid]
        entries_to_delete = delete_feature.to_gff_entry(include_children)
        # delete children
        if include_children and delete_feature.children:
                for child in delete_feature.children[:]:
                    self.pop_feature(child.uid)
        # delete feature from parent
        if delete_feature.puid and delete_feature.puid in self.features:
            self.features[delete_feature.puid].children.remove(delete_feature)
        # delete feature from features
        del self.features[delete_feature.uid]
        # delete feature from lookup
        if delete_feature.aid:
            self.lookup[delete_feature.aid] = [feature for feature in self.lookup[delete_feature.aid] if feature != delete_feature.uid]
            if not self.lookup[delete_feature.aid]:
                del self.lookup[delete_feature.aid]
        return entries_to_delete

    def add_feature(
        self, 
        feature : GFeature, 
        include_children : bool = True
    ) -> str:
        if feature.uid in self.features:
            print("WARNING : duplicate feature already exists and will be overwritten")
        
        self.features[feature.uid] = feature
        if feature.aid:
            if feature.aid in self.lookup:
                self.lookup[feature.aid].append(feature.uid)
            else:
                self.lookup[feature.aid] = [feature.uid]
                
        if feature.paid:
            if feature.paid in self.lookup:
                puid = self.get_uid(feature.paid, feature)
                # feature.set_parent_uid(puid)
                feature.puid = puid

                # TODO: how does this behave?
                if feature not in self.features[puid].children:
                    self.features[puid].add_a_child(feature)
            else:
                print("WARNING: feature has parent attribute, but the parent could not be found in the annotation")
        if include_children:
            for child in feature.children:
                self.add_feature(child)
        return feature.to_gff_entry(include_children)

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

    # TODO: add option to sort
    def to_gff3(self, fp):
        with open(fp, "w") as f:
            f.write("##gff-version 3\n")
            for feature in self.features.values():
                f.write(feature.to_gff_entry(include_children=False))
    def save_as_gix(self, file_path : str):
        with open(file_path, 'wb') as fh:
            pickle.dump(self, fh)
    
    # TODO: fix this
    def clear(self):
        features = dict()
        lookup = dict()

def load_from_gix(file_path : str):
    try:
        with open(file_path, 'rb') as fh:
            res = pickle.load(fh)
    except Exception as e:
        raise RuntimeError(f"error while loading {file_path} : {e}")
    return res