from sage.categories.posets import Posets
from sage.graphs.graph import Graph
from sage.graphs.graph import DiGraph
from sage.combinat.subset import Subsets

def acyclic_orientations(G):

    def reorder_vertices(G):
        n = G.order()
        ko = n
        k = n
        G_copy = G.copy()
        vertex_labels = [None] * n

        while G_copy.size() > 0:
            min_val = float('inf')
            uv = None
            for u, v, _ in G_copy.edges():
                du = G_copy.degree(u)
                dv = G_copy.degree(v)
                val = (du + dv) / (du * dv)
                if val < min_val:
                    min_val = val
                    uv = (u, v)

            if uv:
                u, v = uv
                vertex_labels[u] = ko
                vertex_labels[v] = ko - 1
                G_copy.delete_vertex(u)
                G_copy.delete_vertex(v)
                ko -= 2

            if G_copy.size() == 0:
                break

        for i, label in enumerate(vertex_labels):
            if label is None:
                vertex_labels[i] = ko
                ko -= 1

        return vertex_labels

    def order_edges(G, vertex_labels):
        n = len(vertex_labels)
        m = 1
        edge_labels = {}

        for j in range(2, n + 1):
            for i in range(1, j):
                if G.has_edge(i, j):
                    edge_labels[(i, j)] = m
                    m += 1

        return edge_labels

    def is_upset_of_poset(Poset, subset, keys):
        for (u, v) in subset:
            for (w, x) in keys:
                if (Poset[(u, v), (w, x)] == 1 and (w, x) not in subset):
                    return False
        return True

    def generate_orientations(globO, starting_of_Ek, m, k, keys):
        # Creating a poset
        Poset = {}
        for i in range(starting_of_Ek, m - 1):
            for j in range(starting_of_Ek, m - 1):
                u, v = keys[i]
                w, x = keys[j]
                Poset[(u, v), (w, x)] = 0

        # Create a new graph to determine reachable vertices
        new_G = DiGraph()
        for i in range(0, starting_of_Ek):
            u, v = keys[i]
            if globO[(u, v)] == 1:
                new_G.add_edge(v, u)
            else:
                new_G.add_edge(u, v)

        for i in range(starting_of_Ek, m - 1):
            u, v = keys[i]
            if not new_G.has_vertex(u):
                new_G.add_vertex(u)
            elif not new_G.has_vertex(v):
                new_G.add_vertex(v)

        if (globO[(k-1, k)] == 1):
            new_G.add_edge(k, k - 1)
        else:
            new_G.add_edge(k-1, k)

        for i in range(starting_of_Ek, m - 1):
            for j in range(starting_of_Ek, m - 1):
                u, v = keys[i]
                w, x = keys[j]
                # w should be reachable from u and v should be reachable from x
                if (new_G.shortest_path_length(u, w) < Infinity and new_G.shortest_path_length(x, v) < Infinity):
                    Poset[(u, v), (w, x)] = 1

        # For each subset of the base set of E_k, check if it is an upset or not 
        upsets = []
        for subset in Subsets(keys[starting_of_Ek:m-1]):
            if (is_upset_of_poset(Poset, subset, keys[starting_of_Ek:m-1])):
                upsets.append(list(subset))

        for upset in upsets:
            for i in range(starting_of_Ek, m - 1):
                u, v = keys[i]
                if (u, v) in upset:
                    globO[(u, v)] = 1
                else:
                    globO[(u, v)] = 0

            yield globO.copy()

    def helper(G, globO, m, k):
        keys = list(globO.keys())
        keys = keys[0:m]

        if m <= 0:
            yield {}
            return

        starting_of_Ek = 0
        for (u, v) in keys:
            if u >= k - 1 or v >= k - 1:
                break
            else:
                starting_of_Ek += 1

        # s is the size of E_k
        s = m - 1 - starting_of_Ek

        # Recursively generate acyclic orientations
        orientations_G_small = helper(G, globO, starting_of_Ek, k - 2)

        # For each orientation of G_k-2, yield acyclic orientations
        for alpha in orientations_G_small:
            for (u, v) in alpha:
                globO[(u, v)] = alpha[(u, v)]

            # Orienting H_k as 1
            globO[(k-1, k)] = 1
            yield from generate_orientations(globO, starting_of_Ek, m, k, keys)

            # Orienting H_k as 0
            globO[(k-1, k)] = 0
            yield from generate_orientations(globO, starting_of_Ek, m, k, keys)

    # Reorder vertices based on the logic in reorder_vertices function
    vertex_labels = reorder_vertices(G)
    
    # Create a new graph with updated vertex labels using SageMath
    new_G = Graph()
    for u, v in G.edges(labels=False):  # Assuming the graph edges are unlabelled
        new_G.add_edge(vertex_labels[u], vertex_labels[v])
        
    G = new_G

    # Order the edges based on the logic in order_edges function
    edge_labels = order_edges(G, vertex_labels)

    # Create globO array
    m = len(edge_labels)
    globO = {}
    for (u, v) in edge_labels:
        globO[(u, v)] = 0
        
    k = len(vertex_labels)
    print(vertex_labels)
    print(edge_labels)

    # Call helper function to get acyclic orientations
    orientations = helper(G, globO, m, k)
    count = 0
    for o in orientations:
        print(o)
        count += 1
    print(count)
        
        
# Example usage
# G = graphs.CompleteGraph(5)
G = Graph([(0, 3), (0, 4), (3, 4), (1, 3), (1, 2), (2, 3), (2, 4)])
G.plot()

# Apply the acyclic_orientations function
acyclic_orientations(G)
