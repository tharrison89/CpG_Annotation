"""
Basic utils for annotation
"""

from sklearn.cluster import KMeans

import annotator.cpp_utils as cpp_utils

def compute_gc_content(sequence, kmer_size):
    """
    ComputeGC content with a sliding window across and entire sequence
        Note: hits a cpp util see cpp_utils/gc_content_calculator.cpp for guts

    :param str sequence: Input sequence
    :param int kmer_size: Desired size of sliding window
    :return list(list(float)): A list of features with the only feature being the gc content at each position
    :raises RuntimeError: When the sequence is too small given the kmer size or the kmer size is not valid
    """

    if kmer_size <= 0:
        raise RuntimeError("Please provide a positive number of the kmer size. {} was provided".format(kmer_size))
    if not sequence or kmer_size > len(sequence):
        raise RuntimeError("kmer size of {} is larger than "
                           "the proved sequence of {}".format(kmer_size, sequence))

    return cpp_utils.compute_gc_content(sequence, kmer_size)

def kmeans(input, number_of_clusters):
    """
    Run kmeans on a set of input data

    :param list(list(float)) input: Input where the inner list will contain all features for an element
    :param int number_of_clusters: Desired number of clusters
    :return list(int): Labels for each input
    """
    kmeans = KMeans(n_clusters=number_of_clusters, random_state=0).fit(input)
    labels = kmeans.labels_
    return list(labels)
