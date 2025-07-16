from mjol.gan import *

"""
update_attributes_rule: key -> biotype, value -> list of tuple (attribute name in original, attribute name in other)
ex: 'gene': [('product', 'product'), ('Dbxref', 'dbxref')]
ex: 'default': [('gene', 'gene)] -> 'default' is applied for any feature not included in update attributes rule
"""
def solve_synonym(annotation:GAn, uid:str, other:GAn, other_uid: str,
                  update_attributes_rule: Dict[str, tuple] = {}, 
                  exclude_attributes: List[str] = []):
    
    old_feature = annotation.get_feature(uid)
    new_feature = other.get_feature(other_uid)
    
    raise NotImplementedError