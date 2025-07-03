from pydantic import BaseModel
from typing import Optional, Dict, List, ForwardRef
from intervaltree import IntervalTree

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
    parent : Optional[str] = None
    id : Optional[str] = None
    children : List[GFeatureRef] # type: ignore

    def __init__(self, **data):
        super().__init__(**data)
        self.id = self._infer('id')
        self.parent = self._infer('parent')
    
    def _infer(self, kwd):
        for k, v in self.attributes.items():
            if kwd in k.lower():
                return v
        return None
    
    def __repr__(self) -> str:
        return f"{self.feature_type}:{self._id},{self.strand},{self.start}-{self.end},{self._hash}"

def load_attributes(s: str, kv_sep: str = '=') -> dict:
    return {
        k.strip(): v.strip()
        for x in s.strip().split(';') if x
        for k, v in [x.strip().split(kv_sep, 1)]
    }