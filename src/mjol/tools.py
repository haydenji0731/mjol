from mjol.gan import *
from typing import Tuple

def set_case_insensitive(d, target_key, new_value):
    for k in d:
        if k.lower() == target_key.lower():
            d[k] = new_value
            return
    raise KeyError(f'{target_key} not found')

"""
update_attributes_rule: key -> biotype, value -> list of tuple (attribute name in original, attribute name in other)
ex: 'gene': [('ID, 'ID'), ('product', 'product'), ('Dbxref', 'dbxref')]
ex: 'default': [('gene', 'gene')] -> 'default' is applied for any feature not included in update attributes rule
NOTE: only attributes for gene-level information can be updated (attributes can only be updated using attributes from gene entries)
"""
def solve_synonym(
    original : GAn, 
    uid : str, 
    new : GAn, 
    new_uid : str,
    update_attributes_rule : Dict[str, List[Tuple[str, str]]] = {}, 
    exclude_attributes: List[str] = []
    ) -> tuple:
    
    old_feature = original.get_feature(uid)
    new_feature = new.get_feature(new_uid)
    old_entries = old_feature.to_gff_entry(include_children=True)

    # delete from map
    original.pop_feature(uid, include_children=False)
    
    if old_feature.feature_type in update_attributes_rule:
        for old_attribute, new_attribute in update_attributes_rule[old_feature.feature_type]:
            if (old_attribute in old_feature.attributes) and (new_attribute in new_feature.attributes):
                old_feature.attributes[old_attribute] = new_feature.attributes[new_attribute]
    elif 'default' in update_attributes_rule:
        for old_attribute, new_attribute in update_attributes_rule['default']:
            if (old_attribute in old_feature.attributes) and (new_attribute in new_feature.attributes):
                old_feature.attributes[old_attribute] = new_feature.attributes[new_attribute]
         
    for key in exclude_attributes:
        if key in old_feature.attributes:
                del old_feature.attributes[key]
    # update IDs
    old_feature._populate_gid()
    original.add_feature(old_feature, include_children=False)


    # do the same for children
    for child in old_feature.children[:]:
        # assign parent uid
        set_case_insensitive(child.attributes, child.pak, old_feature.aid)
        child.gid.parent_uid = old_feature.uid
        solve_synonym(original, child.uid, new, new_uid, 
                                update_attributes_rule, exclude_attributes)

    return (old_entries, old_feature.to_gff_entry(include_children=True))