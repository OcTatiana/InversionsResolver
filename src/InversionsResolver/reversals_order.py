from collections import defaultdict
INF = 10 ** 5


def sign(x):
    return (x > 0) - (x < 0)


def reverse(perm, i1, i2):
    result = perm.copy()
    seq = result[i1:i2 + 1]
    reversed_seq = [(-x) for x in reversed(seq)]
    result[i1:i2 + 1] = reversed_seq
    return result


def find_connected_components(index_arrays):
    n = len(index_arrays)
    graph = {i: set() for i in range(n)}

    for i in range(n):
        set_i = set(index_arrays[i])
        for j in range(i + 1, n):
            set_j = set(index_arrays[j])
            if set_i & set_j:
                graph[i].add(j)
                graph[j].add(i)

    visited = set()
    components = []

    for i in range(n):
        if i not in visited:
            component = []
            queue = [i]
            visited.add(i)

            while queue:
                node = queue.pop(0)
                component.append(node)

                for neighbor in graph[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            component.sort()
            components.append([index_arrays[idx] for idx in component])

    return components


def find_non_overlapping_elements(component):
    sets_list = [set(arr) for arr in component]
    non_overlapping = []

    for i, current_set in enumerate(sets_list):
        overlaps_with_any = False

        for j, other_set in enumerate(sets_list):
            if i != j:
                if current_set & other_set and not (current_set.issubset(other_set) or other_set.issubset(current_set)):
                    overlaps_with_any = True
                    break

        if not overlaps_with_any:
            non_overlapping.append(component[i])

    non_overlapping.sort(key=len)

    return non_overlapping


def get_ids_for_perm(components):
    result = []
    for comp_idx, component in enumerate(components):
        comp_id = defaultdict(list)
        if len(component) == 1:
            comp_id[len(component[0])].append(component[0])

        else:
            count = 0
            nested = find_non_overlapping_elements(component)
            for sub in component:
                if sub in nested:
                    comp_id[len(sub)].append(sub)
                else:
                    comp_id[INF + count].append(sub)
                    count += 1

        result.append(comp_id)

    return result

def get_perm_for_image(perm, id_list):
    perm = [perm[0]]
    signs = {abs(i): sign(i) for i in perm[0]}

    components = find_connected_components(id_list)
    result = get_ids_for_perm(components)

    max_iter = max(len(list(comp)) for comp in result)

    for iteration in range(max_iter):
        tmp_perm = perm[-1].copy()
        for comp in result:
            lens = sorted(list(comp))
            if iteration >= len(lens):
                continue
            else:
                length = lens[iteration]
                if length >= INF:
                    id1, id2 = min(tmp_perm.index(i * signs[i]) for i in comp[length][0]), max(
                        tmp_perm.index(i * signs[i]) for i in comp[length][0])
                    tmp_perm = reverse(tmp_perm, id1, id2)
                    for i in comp[length][0]:
                        signs[i] *= -1
                else:
                    for step in comp[length]:
                        id1, id2 = min(tmp_perm.index(i * signs[i]) for i in step), max(
                            tmp_perm.index(i * signs[i]) for i in step)
                        tmp_perm = reverse(tmp_perm, id1, id2)
                        for i in step:
                            signs[i] *= -1
        perm.append(tmp_perm)

    return perm
