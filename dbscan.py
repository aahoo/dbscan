from time import time
from collections import Counter, defaultdict


def dbscan(mat, eps=.15, min_pts=2, report=False):
    '''
    Density-based spatial clustering of applications with noise (DBSCAN).

    :param mat: scipy csr_matrix format
    :param eps: distance epsilon
    :param min_pts: minimum number of points required to form a dense region
    :param report: en/disable measuring duration of each process and reporting
    :return: returns clusters
    '''
    # --------------------------------------------------------------------------

    if report:
        print 'Converting to dictionary',
        init_t = time()

    # @condensed:
    # - keys: frozenset of the index of the nonzero entries of the matrix rows
    # - values: list of the original row indices
    # Note: frozenset is hashable and can be used as a key unlike regular python set.
    # This is useful when obtaining the intersection in Jaccard Distance
    # calculation.
    condensed = defaultdict(list)

    # Note: accessing raw data (indices and indptr attributes of csr_matrix)
    # had 30x performance boost compared to normal indexing/slicing/iteration
    # of sparse matrices.
    j = 0
    for i, k in enumerate(mat.indptr[1:]):
        condensed[frozenset(mat.indices[j:k])].append(i)
        j = k

    if report:
        print 'is done in %d ms.' % ((time() - init_t) * 1000)
        # ----------------------------------------------------------------------
        print 'Constructing adjacency matrix',
        t = time()

    # initializing the adjacency matrix
    adjacency = [[] for i in range(len(condensed))]

    keys = condensed.keys()
    # sorting the points based on number of elements they contains
    sorted_lens = Counter({i: len(key) for i, key in enumerate(keys)}).most_common()

    # constructing adjacency matrix
    for k, (i, i_len) in enumerate(sorted_lens):
        # lengths[k + 1:]: iterating over half-diagonal, and also excluding the
        # diagonal itself since they are alway true
        for j, j_len in sorted_lens[k + 1:]:
            # initial check assumming the best case scenario in which one of
            # the points is subset of the other one. If this neighbor doesn't
            # pass the test, the rest won't either because lengths are sorted
            # so we BREAK the loop. This had 3x performance boost compared to
            # unsorted length case.
            if j_len < (1 - eps) * i_len:
                break

            # calculating proximity based on Jaccard distance
            intersection = len(keys[i] & keys[j])

            # from math: a U b = a + b - a I b
            if 1 - 1.0 * intersection / (i_len + j_len - intersection) <= eps:
                adjacency[i].append(j)
                adjacency[j].append(i)

    if report:
        print 'is done in %d ms.' % ((time() - t) * 1000)
        # ----------------------------------------------------------------------
        print 'DBSCANning',
        t = time()

    # constructing noise cluster and marking them as visited
    visited = [False] * len(condensed)
    clusters = [[]]
    for i, neighbors in enumerate(adjacency):
        if len(neighbors) + len(condensed[keys[i]]) < min_pts:
            clusters[0].append(i)
            visited[i] = True

    # searching for clusters
    for i, neighbors in enumerate(adjacency):
        if not visited[i]:
            visited[i] = True
            clusters.append(set([i]))
            for j in neighbors:
                if not visited[j]:
                    visited[j] = True
                    neighbors.extend([k for k in adjacency[j] if not visited[k]])

                # index -1 is the last item in the list
                clusters[-1].add(j)

    if report:
        print 'is done in %d ms.' % ((time() - t) * 1000)
        # ----------------------------------------------------------------------
        print 'Restoring original point indices',
        t = time()

    # restoring original indices that was lost when condensing the matrix in
    # the beginning
    clusters = [[i for j in c for i in condensed[keys[j]]] for c in clusters]

    if report:
        print 'is done in %d ms.' % ((time() - t) * 1000)
        print '\nProgram finished after %d ms.\n' % ((time() - init_t) * 1000)

        print 'Total number of clusters:', len(clusters)
        print 'Largest cluster size is:', max(map(len, clusters[1:]))
        print 'Noise size:', len(clusters[0])
#         print 'Number of members in each cluster:', map(len, clusters[1:])

    return clusters