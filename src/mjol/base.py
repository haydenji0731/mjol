from pydantic import BaseModel, Field
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
    children : List[GFeatureRef] # type: ignore
    divider : Optional[int] = None
    attributes : Dict[str, str]
    parent : Optional[str] = None
    id : Optional[str] = None

    def __init__(self, **data):
        super().__init__(**data)
        # TODO: fix this
        self.id = self._infer('id')
        self.parent = self._infer('parent')
    
    def _infer(self, kwd):
        for k, v in self.attributes.items():
            if kwd in k.lower():
                return v
        return None
    
    def __repr__(self) -> str:
        return f"{self.feature_type}:{self.id},{self.strand},{self.start}-{self.end}"


    def add_a_child(self, child, divide : bool = True):
        if divide:
            self.divider = len(self.children)
        self.children.append(child)

    def to_gff_entry(self, children=False) -> str:
        attributes_str =";".join(f"{key}={value}" for key, value in self.attributes.items())
        gff_entry = "\t".join([
            str(x) if x is not None else '.' for x in [
                self.chr, self.src, self.feature_type, self.start, self.end, self.score, self.strand, self.frame, attributes_str
            ]
        ])
        if children:
            gff_entry_children = [child.to_gff_entry(children=True) for child in self.children]
            return gff_entry + '\n' + ''.join(gff_entry_children)
        else:
            return gff_entry + '\n'


class IntervalGFeature(GFeature):
    itree: IntervalTree = Field(default_factory=IntervalTree)
    tie_breaking_pad : float = 0.1

    model_config = {
        "arbitrary_types_allowed": True
    }

    def add_a_child(self, child : GFeature):
        self.children.append(child)
        if child.start == child.end:
            self.itree[child.start : child.end + self.tie_breaking_pad] = child
        else:
            self.itree[child.start : child.end] = child
    
    def query_itree(self, start, end):
        if start == end:
            res = self.itree[start : end + self.tie_breaking_pad]
        else:
            res = self.itree[start : end]
        return res

def load_attributes(s: str, kv_sep: str = '=') -> dict:
    return {
        k.strip(): v.strip()
        for x in s.strip().split(';') if x
        for k, v in [x.strip().split(kv_sep, 1)]
    }