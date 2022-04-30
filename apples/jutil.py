from numpy import pad
from apples.support.find_support import get_support, get_support_all

def join_jplace(lst):
    result = lst[0]
    if len(lst) == 1:
        if result["placements"][0]["p"][0][0] == -1:
            result["placements"] = []
    else:
        for i in range(1,len(lst)):
            if lst[i]["placements"][0]["p"][0][0] != -1:
                result["placements"] = result["placements"] + lst[i]["placements"]
    return result

def join_jplace_support(lst):
    support = get_support(lst)
    result = lst[0][0]
    result["placements"][0]["p"][0][2] = support[0]
    if len(lst) == 1:
        if result["placements"][0]["p"][0][0] == -1:
            result["placements"] = []
    else:
        for i in range(1, len(lst)):
            if lst[i][0]["placements"][0]["p"][0][0] != -1:
                lst[i][0]["placements"][0]["p"][0][2] = support[i]
                result["placements"] = result["placements"] + lst[i][0]["placements"]
    return result

def join_jplace_support_all(lst, valids, keep_factor, keep_at_most, prioritize_lse):
    support = get_support_all(lst)
    update_valids(valids, support, keep_factor, keep_at_most, prioritize_lse)
    result = {}
    result["placements"] = []
    for i in range(len(lst)):
        temp = {}
        query_name = lst[i][0]["placements"][0]["n"][0]
        temp["n"] = [query_name]
        temp["p"] = valids[query_name]
        result["placements"].append(temp)
    return result

def update_valids(valids, support, keep_factor, keep_at_most, prioritize_lse):
    for query, placements in valids.items():
        for p in placements:
            if p[0] in support[query]:
                p[2] = support[query][p[0]]
            else:
                p[2] = 0

        if prioritize_lse:
            # only take the placement with minimum LSE
            placements = [min(placements, key=lambda x: x[1])]
        else:
            # sort according to support, then LSE
            placements.sort(key=lambda x: (-x[2], x[1])) # descending by support, ascending by LSE

            # filter out anything with support < keep_factor * highest support
            threshold = keep_factor * placements[0][2]
            placements = list(filter(lambda x: x[2] >= threshold, placements))

            # cut length of list to at most keep_at_most
            if len(placements) > keep_at_most:
                placements = placements[:keep_at_most]

        valids[query] = placements



def extended_newick(tree):
    """Newick printing algorithm is based on treeswift"""

    # if tree.root.edge_length is None:
    #     suffix = ''
    # elif isinstance(tree.root.edge_length, int):
    #     suffix = ':%d' % tree.root.edge_length
    # elif isinstance(tree.root.edge_length, float) and tree.root.edge_length.is_integer():
    #     suffix = ':%d' % int(tree.root.edge_length)
    # else:
    #     suffix = ':%s' % str(tree.root.edge_length)
    suffix = ''
    strng = _nodeprint(tree.root)
    if tree.is_rooted:
        return '[&R] %s%s;' % (strng, suffix)
    else:
        return '%s%s;' % (strng, suffix)


def _nodeprint(root):
    node_to_str = dict()

    for node in root.traverse_postorder():
        if node.is_leaf():
            if node.label is None:
                node_to_str[node] = ''
            else:
                node_to_str[node] = str(node.label)
        else:
            out = ['(']
            for c in node.children:
                out.append(node_to_str[c])
                if c.edge_length is not None:
                    if isinstance(c.edge_length, int):
                        l_str = str(c.edge_length)
                    elif isinstance(c.edge_length, float) and c.edge_length.is_integer():
                        l_str = str(int(c.edge_length))
                    else:
                        l_str = str(c.edge_length)
                    out.append(':%s' % l_str)
                out.append('{%d}' % c.edge_index)
                out.append(',')
                del node_to_str[c]
            out.pop()  # trailing comma
            out.append(')')
            if node.label is not None:
                out.append(str(node.label))
            node_to_str[node] = ''.join(out)
    return node_to_str[root]