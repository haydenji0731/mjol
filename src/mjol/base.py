from pydantic import BaseModel
from typing import Optional, Dict, List, ForwardRef

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
    children : List[GFeatureRef] # type: ignore
    attributes : Dict[str, str]
    parent : Optional[str] = None
    id : Optional[str] = None

    def __init__(self, **data):
        super().__init__(**data)
        self.id = self._infer('id')
        self.parent = self._infer('parent')
    
    def _infer(self, kwd):
        for k, v in self.attributes.items():
            if kwd in k.lower():
                return v
        return None

    def build_itree_from_children(self):
        raise NotImplementedError
    
    def __repr__(self) -> str:
        return f"{self.feature_type}:{self.id},{self.strand},{self.start}-{self.end}"

def load_attributes(s: str, kv_sep: str = '=') -> dict:
    return {
        k.strip(): v.strip()
        for x in s.strip().split(';') if x
        for k, v in [x.strip().split(kv_sep, 1)]
    }