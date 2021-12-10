from collections import deque
import numpy as np

class Node:
    def __init__(self, name='', edge_len = 0, parent = None):
        self.name = name
        self.length = edge_len
        self.parent = parent
        self.children = []
        self.SNVs = []
        self.CNVs = []
        self.ancestral_SNVs = []
        self.ancestral_CNVs = []
        self.genome = None

    def __str__(self):
        if self.name:
            return str(self.name)
        else:
            return ''
    
    def write_newick(self, terminate=True, format=0):
        if self.is_leaf():
            if format == 0:
                return self.name
            elif format == 1:
                return self.name + ':' + str(self.length)
        else:
            newick_str = '('
            for child in self.children:
                newick_str += child.write_newick(terminate=False, format=format) + ','
            newick_str = newick_str[:-1]
            newick_str += ')'
            if format == 0 or terminate:
                newick_str += self.name
            elif format == 1:
                newick_str += self.name + ':' + str(self.length)
        if terminate:
            return newick_str + ';'
        else:
            return newick_str

    def set_name(self, name):
        self.name = name

    def set_len(self, edge_len):
        self.length = float(edge_len)

    def set_parent(self, parent):
        self.parent = parent

    def set_child(self, child):
        self.children.append(child)
        child.parent = self

    def set_sibling(self, sibling):
        self.sibling.append(sibling)

    def is_leaf(self):
        return len(self.children) == 0
    
    def is_root(self):
        return self.parent is None

    def get_root(self):
        root = self
        while root.parent is not None:
            root = root.parent
        return root

    def add_child(self, name='', edge_len = 0):
        child = Node(name = name, edge_len = float(edge_len), parent = self)
        self.children.append(child)
        return child

    def iter_postorder(self):
        visit_queue = deque()
        return_queue = deque()
        visit_queue.append(self)

        while visit_queue:
            node = visit_queue.pop()
            return_queue.append(node)
            if not node.is_leaf():
                visit_queue.extend(node.children)
        
        while return_queue:
            node = return_queue.pop()
            yield node

    #level order traversal, i.e. breadth first search from the root. Does include the root itself.
    def iter_descendants(self):
        nodes = deque()
        nodes.append(self)

        while nodes:
            node = nodes.popleft()
            nodes.extend(node.children)
            #if node != self:
            yield node

    def get_leaves(self):
        return [n for n in self.iter_postorder() if n.is_leaf()]

    def get_height(self):
        if self.is_leaf():
            return 0
        else:
            return max([1 + child.get_height() for child in self.children])

    def get_total_branchlen(self):
        return sum([node.length for node in self.iter_descendants()])

class Tree:
    def __init__(self, root=None, newick=None):
        self.root = root
        if newick:
            self.root = str_to_newick(newick)
        self.cell_names = [node.name for node in self.iter_leaves()]

    def print_newick(self, format=0):
        newick_str = self.root.write_newick(format=format)
        print(newick_str)

    def iter_postorder(self):
        return self.root.iter_postorder()
    
    def iter_descendants(self):
        return self.root.iter_descendants()
    
    def iter_leaves(self):
        for node in self.iter_descendants():
            if node.is_leaf():
                yield node

    def set_leaf_names(self):
        count = 1
        for leaf in self.iter_leaves():
            leaf.name = 'cell' + str(count)
            count += 1
        self.cell_names = [node.name for node in self.iter_leaves()]

    def get_tree_height(self):
        return self.root.get_height()
    
    def get_total_branchlen(self):
        return self.root.get_total_branchlen()

def str_to_newick(newick_str):
    split_str = newick_str[:-1].split(',')

    cur_node = None
    nodes = []

    for chunk in split_str:
        while chunk[0] == '(':
            new_node = Node()
            if cur_node:
                cur_node.set_child(new_node)
            cur_node = new_node
            chunk = chunk[1:]
        rest = chunk.split(')')
        if ':' in rest[0]:
            idx = rest[0].index(':')
            cur_node.add_child(name=rest[0][:idx], edge_len=rest[0][idx+1:])
        else:
            cur_node.add_child(name=rest[0])
        if len(rest) > 1:
            for part in rest[1:]:
                if ':' in part:
                    idx = part.index(':')
                    cur_node.set_name(part[:idx])
                    cur_node.set_len(part[idx+1:])
                else:
                    cur_node.set_name(part)
                    cur_node.set_len(0)
                if not cur_node.is_root():
                    cur_node = cur_node.parent
    
    return cur_node

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