import heapq


def overlap(r1, r2):
    i1, j1 = r1
    i2, j2 = r2
    return not (j1 < i2 or j2 < i1)


def build_dependency_graph(reversals):
    n = len(reversals)
    graph = {i: set() for i in range(n)}

    for i in range(n):
        for j in range(i+1, n):
            if overlap(reversals[i], reversals[j]):
                graph[i].add(j)

    return graph


def canonical_order(reversals, graph):
    n = len(reversals)
    indegree = [0] * n

    for u in graph:
        for v in graph[u]:
            indegree[v] += 1

    heap = []
    for i in range(n):
        if indegree[i] == 0:
            heapq.heappush(heap, (reversals[i], i))

    order = []
    while heap:
        _, u = heapq.heappop(heap)
        order.append(u)

        for v in graph[u]:
            indegree[v] -= 1
            if indegree[v] == 0:
                heapq.heappush(heap, (reversals[v], v))

    return [reversals[i] for i in order]


def canonicalize(reversals):
    graph = build_dependency_graph(reversals)
    return canonical_order(reversals, graph)
