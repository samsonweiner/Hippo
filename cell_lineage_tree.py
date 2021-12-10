#Tree with random topology by attaching two random leaves. Branch lengths are left undefined.
def gen_random_topology(m, names = None):
    root = Node()
    leaves = [root]

    while len(leaves) < m:
        idx = np.random.randint(0, len(leaves))
        node = leaves.pop(idx)
        c1 = node.add_child()
        c2 = node.add_child()
        leaves.append(c1)
        leaves.append(c2)
    
    count = 0
    for n in root.get_leaves():
        if names:
            n.name = names[count]
        else:
            n.name = 'leaf' + str(count)
        count += 1

    t = Tree(root=root)
    return t

# Creates branch lengths so that each leaf is equally distant from the root.
def add_branchlen_ultrametric(tree):
    node_max_depths = {}
    for node in tree.iter_postorder():
        if node.is_leaf():
            node_max_depths[node] = 1
        else:
            node_max_depths[node] = max([node_max_depths[c] for c in node.children]) + 1
    
    tree_length = float(tree.get_tree_height())
    node_dists = {tree.root: 0.0}
    for node in tree.iter_descendants():
        node.length = (tree_length - node_dists[node.parent]) / node_max_depths[node]
        node_dists[node] = node.length + node_dists[node.parent]
    
# Given a tree topology with existing branch lengths, computes branch lengths uniformly with a random deviation proportional to the tree height.
def add_branchlen_deviation(tree, shape=1):
    for node in tree.iter_descendants():
        x = np.random.gamma(shape)
        node.length = node.length * x