from pydantic import BaseModel, Field
from typing import Optional, Dict, List, ForwardRef
import hashlib
import pickle

class GId(BaseModel):
    uid : Optional[str] = None
    aid : Optional[str] = None
    parent : Optional[str] = None
    parent_uid : Optional[str] = None

GFeatureRef = ForwardRef("GFeature")

class GFeature(BaseModel):
    chr : str
    src : str
    feature_type : str
    start : int
    end : int
    score : Optional[float]
    strand : str # ['.', '-', '+']
    frame : str # ['.', '0', '1', '2']
    attributes : Dict[str, str]
    children : List[GFeatureRef] = Field(default_factory=list) # type: ignore
    iak : str
    pak : str
    gid : GId = Field(default_factory=GId)

    def __init__(self, **data):
        super().__init__(**data)
        self._populate_gid()
    
    def __repr__(self) -> str:
        return f"{self.feature_type}:{self.aid or ''}:{self.uid},{self.chr},{self.strand},{self.start}-{self.end}"
    
    def _populate_gid(self):
        self.gid.uid = self._assign_uid()
        self.gid.aid = self._infer(self.iak)
        self.gid.parent = self._infer(self.pak)

    def _infer(self, ak):
        for k, v in self.attributes.items():
            if ak.lower() == k.lower():
                return v
        return None

    def _assign_uid(self):
        s = f'{self.chr}.{self.start}.{self.end}.{self.strand}.{self.feature_type}'
        uid = hashlib.sha256(s.encode()).hexdigest()
        return uid
    
    def add_a_child(self, child):
        self.children.append(child)
    
    def set_parent_uid(self, puid):
        self.gid.parent_uid = puid
    
    def to_gff_entry(self, include_children : bool = False):
        attributes_str =";".join(f"{key}={value}" for key, value in self.attributes.items())
        entry = "\t".join([
            str(x) if x is not None else '.' for x in [
                self.chr, self.src, self.feature_type, self.start, self.end, self.score, self.strand, self.frame, attributes_str
            ]
        ])
        if include_children:
            children_entry = [child.to_gff_entry(include_children=True) for child in self.children]
            return entry + '\n' + ''.join(children_entry)
        return entry + '\n'
    
    def calc_sim(self, other):
        score = 0
        if self.chr == other.chr:
            score += 1
        if self.strand == other.strand:
            score += 1
        score += abs(self.start - other.start)
        score += abs(self.end - other.end)
        return score

    @property
    def uid(self):
        return self.gid.uid
    
    @property
    def aid(self):
        return self.gid.aid
    
    @property
    def parent(self):
        return self.gid.parent

    @property
    def parent_uid(self):
        return self.gid.parent_uid