# dist_utils.py

def jaccard(a, b):
    """
    Compute the Jaccard distance between two sets. Note that Jacc({}, {}) = 0.
    :param a: a set
    :param b: another set
    :return: Jacc(a, b)
    """
    union = len(a | b)
    intersection = len(a & b)
    return (1 - intersection / union) if union > 0 else 0


def ancestor_sets(newick):
    """
    Find the ancestor sets of every mutation in tree represented in Newick format.
    :param newick: a Newick string representation of a tree
    :return: a dict of each mutation's set of ancestors
    """
    ancestors = dict()
    _ancestor_sets(newick.replace(' ', ''), ancestors, set())
    return ancestors


def _ancestor_sets(newick, ancestors, current_ancestors):
    """
    Recursive helper function for ancestor_sets().
    :param newick: a Newick string representation of a subtree
    :param ancestors: a pointer to the dict being computed
    :param current_ancestors: the mutations ancestral to every mutation in the current subtree
    """
    last_paren_index = newick.rfind(')')
    label_string = newick[last_paren_index + 1:]

    if label_string[0] == '{':
        labels = set(label_string[1:-1].split(','))
    else:
        labels = {label_string}

    for label in labels:
        ancestors[label] = labels | current_ancestors

    if last_paren_index != -1:
        for subtree_newick in outer_comma_split(newick[1:last_paren_index]):
            _ancestor_sets(subtree_newick, ancestors, labels | current_ancestors)


def outer_comma_split(newick):
    """
    Split a Newick subtring on commas, ignoring those contained in parentheses and brackets.
    :param newick: the string to split
    """
    chunk_start = 0
    parens = 0
    brackets = 0

    for i, char in enumerate(newick):
        if char == '(':
            parens += 1
        elif char == ')':
            parens -= 1
        elif char == '{':
            brackets += 1
        elif char == '}':
            brackets -= 1
        elif char == ',' and parens == 0 and brackets == 0:
            yield newick[chunk_start:i]
            chunk_start = i + 1

    yield newick[chunk_start:]