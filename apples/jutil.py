from numpy import pad
from apples.support.find_support import get_support, get_support_all

def join_jplace(lst):
    """
    Join the list of jplace files into a single jplace file.

    Parameters:
    lst (list): A list of jplace files.

    Returns:
    dict: The joined jplace dictionary.
    """
    result = lst[0]
    if len(lst) == 1:
        if result['placements'][0]['p'][0][0] == -1:
            result['placements'] = []
    else:
        for i in range(1, len(lst)):
            if lst[i]['placements'][0]['p'][0][0] != -1:
                result['placements'] = result['placements'] + lst[i]['placements']
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
                result["placements"] += lst[i][0]["placements"]
    return result

def join_jplace_support_all(lst, valids, keep_factor, keep_at_most, prioritize_lse):
    """
    Join the list of jplace files into a single jplace file. Include support values of all placements too.

    Parameters:
    lst (list): A list of jplace files.

    Returns:
    dict: The joined jplace dictionary.
    """
    support = get_support_all(lst)
    update_valids(valids, support, keep_factor, keep_at_most, prioritize_lse)
    result = {}
    result["placements"] = []
    for query in lst:
        temp = {}
        temp["n"] = [query]
        temp["p"] = valids[query]

        # only append if there is any valid placement
        if len(temp["p"]) > 0:
            result["placements"].append(temp)
    return result

def update_valids(valids, support, keep_factor, keep_at_most, prioritize_lse):
    """
    Comment required
    """
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

        # cut out -1 branches
        valids[query] = list(filter(lambda x: x[0] != -1, valids[query]))

def extended_newick(tree):
    """
    Generate an extended Newick string representation of a tree using treeswift.

    This function prepares a suffix as an empty string and uses the recursive helper function _nodeprint
    to create a string representation of the tree with edge indices. It returns the Newick string
    along with additional formatting to denote if the tree is rooted.

    Args:
        tree: A Tree object to be converted into an extended Newick format.

    Returns:
        str: An extended Newick string representation of the tree.
"""

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
    """
    Generate a string representation of a tree from its root node.
    This code defines a function _nodeprint that takes a root node and
    returns a dictionary node_to_str containing string representations of
    nodes in a tree. It uses postorder traversal to process each node and its
    children, and constructs a string representation for each node based on
    its label, edge length, and **edge index**.

    Args:
        root: The root node of the tree.

    Returns:
        A string representation of the tree.
    """
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
