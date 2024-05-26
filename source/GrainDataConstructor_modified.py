from collections import defaultdict
from operator import itemgetter

def group_topology(topology):
    """Groups all atoms that belong to a certain GB, TL, QP or IA together without duplication."""
    new_topology = defaultdict(list)
    for parent_grains, atom_position in topology:
        key = tuple(sorted(parent_grains))
        new_topology[key].append(atom_position)
    return sorted(new_topology.items(), key=lambda x: int(x[0][0]))

def calculate_grain_number(grain_id_list):
    """Computes the number of Grains that are in our microstructure."""
    return int(max(grain_id_list, key=itemgetter(1))[1]) + 1

def construct_grain_info(grain_list, X):
    """Creates structured data for each kind of network entity."""
    entities = defaultdict(list)
    for grain in grain_list[:-1]:
        atom_id,modified_grain_id, topology_list= grain
        atom_position = X[atom_id].tolist()
        if modified_grain_id>=0:
            entity_index=0
        elif modified_grain_id<0:
            entity_index=int(abs(modified_grain_id))
        
        entity_type = ['inner_atom', 'grain_boundary', 'triple_line', 'quadraple_point'][entity_index]
        entities[entity_type].append([list(map(int, topology_list)), atom_position])
    return (entities['inner_atom'], entities['grain_boundary'], entities['triple_line'], entities['quadraple_point'])

def order_topologies(modified_grain_id_list, X):
    """Orders the structured network entity data according to the Id of the first grain."""
    inner_atom, grain_boundary, triple_line, quadraple_point= construct_grain_info(modified_grain_id_list, X)
    return {entity: group_topology(data) for entity, data in 
            [('inner_atom', inner_atom), 
             ('grain_boundary', grain_boundary), 
             ('triple_line', triple_line), 
             ('quadraple_point', quadraple_point)]}
