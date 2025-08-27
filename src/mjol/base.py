from pydantic import BaseModel, Field
from typing import Optional, Dict, List, ForwardRef
import hashlib

class GId(BaseModel):
    uid : Optional[str] = None
    aid : Optional[str] = None
    paid : Optional[str] = None # parent_aid
    puid : Optional[str] = None # parent_uid

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
        self.gid.paid = self._infer(self.pak)

    def _infer(self, ak):
        for k, v in self.attributes.items():
            if ak.lower() == k.lower():
                return v
        return None

    def _assign_uid(self):
        """
        hashes the entire gtf/gff line to obtain a uid
        """
        s = f'{self.chr}\t{self.src}\t{self.feature_type}\t{self.start}\t{self.end}\t'
        s += f'{'.' if not self.score else self.score}\t{self.strand}\t{self.frame}\t'
        s += ';'.join([f'{k}={v}' for k, v in self.attributes.items()])
        uid = hashlib.sha256(s.encode()).hexdigest()
        return uid
    
    def add_a_child(self, child):
        self.children.append(child)

    # TODO: remove this func usage everywhere
    def set_parent_uid(self, puid):
        self.gid.puid = puid
    
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
    
    def get_chain(self, feature_type, sort=True):
        chain = []
        for child in self.children:
            if child.feature_type == feature_type:
                chain.append((child.start, child.end))
        if sort:
            return sorted(chain)
        return chain

    # getter, setter methods
    
    @property
    def uid(self):
        return self.gid.uid
    
    @uid.setter
    def uid(self, new_uid):
        self.gid.uid = new_uid
        
    @property
    def aid(self):
        return self.gid.aid

    @aid.setter
    def aid(self, new_aid):
        self.gid.aid = new_aid
    
    @property
    def paid(self):
        return self.gid.paid

    @paid.setter
    def paid(self, new_paid):
        self.gid.paid = new_paid

    @property
    def puid(self):
        return self.gid.puid

    @puid.setter
    def puid(self, new_puid):
        self.gid.puid = new_puid