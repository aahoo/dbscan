
from time import gmtime, strftime, time
from collections import Counter


def dbscan(sparse_mat, eps=.15, min_pts=2, report=False):
    """
    Density-based spatial clustering of applications with noise (DBSCAN).

    :param sparse_mat: Scipy sparse matrix format
    :param eps: distance epsilon
    :param min_pts: minimum number of points required to form a dense region
    :param report: en/disable measuring duration of each process and reporting
    :return: returns nothing
    """
    # --------------------------------------------------------------------------
    if report:
        print "Converting to dictionary",
        init_time = time()

    # converting to lil matrix and getting access to RAW data had 20x
    # performance boost
    data = sparse_mat.tolil()

    # keys are frozenset of the index of the nonzero entries of the matrix
    # rows. frozenset is hashable and can be used as the key unlike regular
    # sets. values are the row index
    dictionary = {}
    for i, row in enumerate(data.rows):
        key = frozenset(row)
        if key in dictionary:
            dictionary[key].add(i)
        else:
            dictionary[key] = set([i])

    # this is used later for optimizing construction of neighbors matrix
    lengths = Counter()
    keys = dictionary.keys()
    for i, key in enumerate(keys):
        lengths[i] = len(key)

    if report:
        print "is done in %.1fs." % (time() - init_time)
        # ----------------------------------------------------------------------
        print "Constructing neighbors matrix",
        t = time()

    # initializing neighbors
    neighbors = [[] for i in range(len(dictionary))]

    # constructing neighbors matrix
    sorted_lengths = lengths.most_common()
    for idx, (i, i_len) in enumerate(sorted_lengths):
        # lengths[idx+1:]: iterating over half-diagonal, and also excluding the
        # diagonal itself since they are alway true
        for i_neibor, i_neibor_len in sorted_lengths[idx + 1:]:
            # initial check assumming the best case scenario in which one of
            # the points is subset of the other one. If this neighbor doesn't
            # pass the test, the rest won't either because lengths are sorted
            # so we BREAK the loop. This had 3x performance boost compared to
            # unsorted length case.
            if i_neibor_len < (1 - eps) * i_len:
                break

            # calculating proximity based on Jaccard distance
            intersection = len(keys[i] & keys[i_neibor])
            # from math: len(a U b) = len(a) + len(b) - len(a I b)
            if 1 - 1.0 * intersection / (i_len + i_neibor_len - intersection) <= eps:
                neighbors[i].append(i_neibor)
                neighbors[i_neibor].append(i)

    if report:
        print "is done in %.1fs." % (time() - t)
        # ----------------------------------------------------------------------
        print "DBSCANning",
        t = time()

    # constructing noise cluster
    visited_noise = [False] * len(dictionary)
    noise = set([])
    for i, i_neighbors in enumerate(neighbors):
        if len(i_neighbors) + len(dictionary[keys[i]]) < min_pts:
            # Restoring original indices of the nodes
            noise.update(dictionary[keys[i]])
            visited_noise[i] = True

    # searching for clusters
    clusters = []
    for i, i_neighbors in enumerate(neighbors):
        if not visited_noise[i]:
            visited_noise[i] = True
            # dictionary[keys[i]] restores the original indices of the nodes
            clusters.append(dictionary[keys[i]])
            for i_neibor in i_neighbors:
                if not visited_noise[i_neibor]:
                    visited_noise[i_neibor] = True
                    for i_neibor_neibor in neighbors[i_neibor]:
                        if not visited_noise[i_neibor_neibor]:
                            i_neighbors.append(i_neibor_neibor)
                # index -1 is the last item in the list
                clusters[-1].update(dictionary[keys[i_neibor]])

    if report:
        print "is done in %.1fs." % (time() - t)
        print "Program finished after %.1fs.\n" % (time() - init_time)

        print "Largest cluster size is:\t", max(map(len, clusters))
        print "Total number of clusters:\t", len(clusters) + 1
        print "Number of members in noise:\t", len(noise)
        # print "Number of members in each cluster:\t", map(len, clusters)
