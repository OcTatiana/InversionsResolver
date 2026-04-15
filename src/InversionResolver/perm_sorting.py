import networkx as nx
from itertools import combinations
import random
import matplotlib.pyplot as plt


def build_points_order(perm):
    # [(-0, +0), (-3, +3), ...]
    pts = []
    for x in perm:
        v = abs(x)
        if x >= 0:
            pts.append(('-', v)); pts.append(('+', v))
        else:
            pts.append(('+', v)); pts.append(('-', v))
    return pts


def find_point_positions(pts):
    pos = {}
    for i, p in enumerate(pts):
        pos.setdefault(p, []).append(i)
    return pos


def arc_endpoints_positions(perm):
    pts = build_points_order(perm)
    posmap = find_point_positions(pts)
    # print(posmap)
    n = len(perm) - 2
    endpoints = []
    for i in range(0, n + 1):
        p_plus_i = ('+', i)
        p_minus_i1 = ('-', i + 1)

        if p_plus_i not in posmap or p_minus_i1 not in posmap:
            raise ValueError(f"Missing point {p_plus_i} or {p_minus_i1} for arc {i}.")

        pa = posmap[p_plus_i][0]
        pb = posmap[p_minus_i1][0]
        endpoints.append((min(pa, pb), max(pa, pb)))
    return endpoints


def build_overlap_graph(perm):
    n = len(perm) - 2
    G = nx.Graph()
    endpoints = arc_endpoints_positions(perm)
    # print(endpoints)
    perm_sorted_abs = sorted(perm, key=abs)
    # print(perm_sorted_abs)
    for i in range(n + 1):
        G.add_node(i)
        G.nodes[i]['oriented'] = ((perm_sorted_abs[i]+0.1) * perm_sorted_abs[i + 1] < 0)
        G.nodes[i]['arc_endpoints'] = endpoints[i]
    for i, j in combinations(range(n + 1), 2):
        a1, a2 = endpoints[i]; b1, b2 = endpoints[j]
        if (a1 < b1 < a2 < b2) or (b1 < a1 < b2 < a2):
            G.add_edge(i, j)
    return G


def draw_overlap_graph(G: nx.Graph, title: str = "Overlap graph"):
    if len(G) == 0:
        print("[draw] empty graph")
        return
    pos = nx.circular_layout(G) if len(G) > 1 else {list(G.nodes())[0]: (0, 0)}
    node_colors = ['black' if G.nodes[v]['oriented'] else 'white' for v in G.nodes()]
    nx.draw(G, pos, with_labels=False, node_color=node_colors, edgecolors='black', node_size=700)

    labels = {v: f"v{v}" for v in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_color='red')
    plt.title(title)
    plt.show()


def has_unoriented_component(G: nx.Graph):
    for comp in nx.connected_components(G):
        if len(comp) <= 1:
            continue
        if not any(G.nodes[v]['oriented'] for v in comp):
            return True
    return False


def arc_to_reversal_interval(perm, arc_index):
    pts = build_points_order(perm)
    posmap = find_point_positions(pts)
    endpoints = arc_endpoints_positions(perm)
    a, b = endpoints[arc_index]
    n = len(perm)

    labels_with_point_inside = set()
    for k in range(0, n):
        p_plus = ('+', k)
        p_minus = ('-', k)

        inside = False
        if p_plus in posmap:
            for pa in posmap[p_plus]:
                if a < pa < b:
                    inside = True
                    break

        if (not inside) and (p_minus in posmap):
            for pb in posmap[p_minus]:
                if a < pb < b:
                    inside = True
                    break
        if inside:
            labels_with_point_inside.add(k)

    labels_with_point_inside.discard(0)
    labels_with_point_inside.discard(n-1)

    if not labels_with_point_inside:
        return None

    positions = [i for i, v in enumerate(perm) if abs(v) in labels_with_point_inside]
    if not positions:
        return None

    L = min(positions)
    R = max(positions)
    return (L, R)


def apply_reversal(perm, L, R):
    new = perm[:L] + [-x for x in perm[L:R+1][::-1]] + perm[R+1:]
    return new


def local_complement(G: nx.Graph, v):
    if v not in G:
        raise KeyError(f"Vertex {v} not in graph")

    H = G.copy()
    neighbors = list(H.neighbors(v)) + [v]

    for u1, u2 in combinations(neighbors, 2):
        if H.has_edge(u1, u2):
            H.remove_edge(u1, u2)
        else:
            H.add_edge(u1, u2)

    for u in neighbors:
        H.nodes[u]['oriented'] = not H.nodes[u]['oriented']

    return H


def is_isolated_in_subgraph(G, v, Q):
    for nbr in G.neighbors(v):
        if nbr in Q:
            return False
    return True


def has_good(G, Q):
    for u in Q:
        if G.nodes[u].get("oriented") and not is_isolated_in_subgraph(G, u, Q):
            return True
    return False


def get_any_good(G, Q):
    for u in Q:
        if G.nodes[u].get("oriented") and not u in nx.isolates(G):
            return u
    return None


def do_good(G, Q):
    S = []

    while True:
        v = get_any_good(G, Q)
        if v is None:
            break
        G = local_complement(G, v)
        S.append(v)
        Q = [x for x in Q if not x in nx.isolates(G)]
        random.shuffle(Q)

    return S, Q


def recover(G, Sf, Q):
    S1 = list(Sf)
    S2 = []
    tmpG = [G]

    for v in S1:
        tmpG.append(local_complement(tmpG[-1], v))

    while S1 and not has_good(tmpG[-1], Q):
        v = S1.pop()
        tmpG.pop()
        S2.insert(0, v)
        # Q.add(v)
    return S1, S2, tmpG[-1], Q


def sort_graph(G, Q):
    Sf, Q = do_good(G, Q)
    # print("Made good: ", Sf, Q)
    Sb = []

    while Sf:
        S1, S2, tmpG, Q = recover(G, Sf, Q)
        # print("After recover: ", S1, S2, Q)
        Sprime, Q = sort_graph(tmpG, Q)
        # print("S': ", Sprime, Q)

        if len(Sprime) % 2 == 0:
            Sb = Sprime + S2 + Sb

        else:
            if S2:
                S2_without_first = S2[1:]
            else:
                S2_without_first = []
            Sb = Sprime + S2_without_first + Sb

        Sf = S1
        # print("end: ", Sf, Sb)

    return Sf + Sb, Q


def init_perm(perm):
    n = len(perm) - 2
    G = nx.Graph()
    endpoints = arc_endpoints_positions(perm)
    # print(endpoints)
    perm_sorted_abs = sorted(perm, key=abs)
    # print(perm_sorted_abs)
    for i in range(n + 1):
        G.add_node(i)
        G.nodes[i]['oriented'] = ((perm_sorted_abs[i]+0.1) * perm_sorted_abs[i + 1] < 0)
        G.nodes[i]['arc_endpoints'] = endpoints[i]
    for i, j in combinations(range(n + 1), 2):
        a1, a2 = endpoints[i]; b1, b2 = endpoints[j]
        if (a1 < b1 < a2 < b2) or (b1 < a1 < b2 < a2):
            G.add_edge(i, j)
    return G


def sort_perm(G, Q):
    Sf, Q = do_good(G, set(Q))
    # print("Made good: ", Sf, Q)
    Sb = []

    while Sf:
        S1, S2, tmpG, Q = recover(G, Sf, Q)
        # print("After recover: ", S1, S2, Q)
        Sprime, Q = sort_graph(tmpG, Q)
        # print("S': ", Sprime, Q)

        if len(Sprime) % 2 == 0:
            Sb = Sprime + S2 + Sb

        else:
            if S2:
                S2_without_first = S2[1:]
            else:
                S2_without_first = []
            Sb = Sprime + S2_without_first + Sb

        Sf = S1
        # print("end: ", Sf, Sb)

    return Sf + Sb, Q


def find_reversal_boundaries(p1, p2):
    boundaries = []
    if p1 == p2:
        return boundaries

    diff = [i for i in range(len(p1)) if p1[i] != p2[i]]
    if not diff:
        return boundaries

    left, right = diff[0], diff[-1]
    boundaries.append((abs(p1[left]), abs(p1[right])))
    return boundaries


def get_reversal_sequences(perm, N = 10000):
    sorted_perm = sorted([abs(x) for x in perm])

    change_sign = False
    if all(x >= 0 for x in perm):
        change_sign = True

    sequences = []

    for _ in range(N):
        if change_sign:
            el_id = random.randint(1, len(perm)-2)
            perm[el_id] = -perm[el_id]

            G = build_overlap_graph(perm)
            perm[el_id] = -perm[el_id]
            Q = [i for i in range(len(perm)-1)]
            random.shuffle(Q)
            seq, final_Q = sort_graph(G, Q)
            if seq not in sequences:
                sequences.append((seq, el_id))

        else:
            G = build_overlap_graph(perm)
            Q = [i for i in range(len(perm)-1)]
            random.shuffle(Q)
            seq, final_Q = sort_graph(G, Q)
            if seq not in sequences:
                sequences.append((seq, -1))

    perm_global_list = []
    id_global_list = []
    block_global_list = []
    min_len = min(len(x[0]) for x in sequences)

    for seq, el_id in sequences:
        if perm_global_list and len(seq) > len(perm_global_list[0]) - 1:
            continue

        id_list = []
        block_list = []

        perm_list = [perm]
        if change_sign:
            perm[el_id] = -perm[el_id]
            perm_list.append(perm.copy())

        tmp_perm = perm.copy()

        if change_sign:
            perm[el_id] = -perm[el_id]

        for i in seq:
            try:
                L, R = arc_to_reversal_interval(tmp_perm, i)
                id_list.append((L, R))
                block_list.append(sorted([abs(x) for x in tmp_perm[L:R + 1]]))
                tmp_perm = apply_reversal(tmp_perm, L, R)
                perm_list.append(tmp_perm)

            except TypeError:
                break

        block_list.sort()
        # print(block_list)
        if perm_list not in perm_global_list and block_list not in block_global_list and perm_list[-1] == sorted_perm:
            perm_global_list.append(perm_list)
            id_global_list.append(id_list)
            block_global_list.append(block_list)
        # order = canonicalize(id_list)
        # if perm_list not in perm_global_list and order not in id_global_list and perm_list[-1] == sorted_perm:
        #     perm_global_list.append(perm_list)
        #     id_global_list.append(order)
    # print(*block_global_list, sep = "\n")

    return perm_global_list, block_global_list