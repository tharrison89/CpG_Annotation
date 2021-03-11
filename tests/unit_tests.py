"""
Unit tests for this project
"""

import unittest

import annotator.annotation_utils as annotation_utils
from annotator.annotator_base import StateInterval
from annotator.kmeans_gc_content_annotator import KmeansGcContentAnnotator

class GCKmerizerTests(unittest.TestCase):
    """
    Tests for the GC kmerizer
    """

    def testHappyPath(self):
        """
        Test standard path
        """
        INPUT_STR = "AGGCCTA"
        INPUT_KMER_LEN = 2
        EXPECTED_RESULTS = [[0.5], [1.0], [1.0], [1.0], [0.5], [0.0]]
        results = annotation_utils.compute_gc_content(INPUT_STR, INPUT_KMER_LEN)
        self.assertListEqual(results, EXPECTED_RESULTS)

    def testInputShorterThanKmer(self):
        """
        Check for the case where the input sequece is smaller than the kmer size
        """
        INPUT_STR = "AGGCCTA"
        INPUT_KMER_LEN = 200
        with self.assertRaises(RuntimeError):
            annotation_utils.compute_gc_content(INPUT_STR, INPUT_KMER_LEN)

    def testInvalidKmer(self):
        """
        Test for an invalid kmer length
        """
        INPUT_STR = "AGGCCTA"
        INPUT_KMER_LEN = -2
        with self.assertRaises(RuntimeError):
            annotation_utils.compute_gc_content(INPUT_STR, INPUT_KMER_LEN)

class KmeansTests(unittest.TestCase):
    """
    Test kmeans calculation
    """

    def testOneDimensionsKmeans(self):
        """
        Test one dimensional kmeans
        """
        INPUT = [[0.01], [0.02], [1.0], [1.1], [1.3], [0.03]]
        NUMBER_OF_CLUSTERS = 2
        results = annotation_utils.kmeans(INPUT, NUMBER_OF_CLUSTERS)
        EXPECTED_LABELS = [1, 1, 0, 0, 0, 1]
        self.assertListEqual(results, EXPECTED_LABELS)

class GcContentKmeansTest(unittest.TestCase):
    """
    Tests for the GC content anotator
    """

    def intervalToTupple(self, interval):
        """
        Helper function to represent the interval object as a tuple for easier comparison

        :param GCScoredInterval interval: Input interval
        :return object, int, int: The label, start, and stop respectively.
        """
        return (interval.label, interval.start, interval.stop)

    def testEasyExampleNoSmoothing(self):
        """
        Test the states without smoothing
        """

        INPUT_SEQUENCE = "GCGCCCCGCAGCGCGATATATATATATAATATGCATATATATATATGCGCGCGCGCGGCGCGCGCGC"

        EXPECTED_RESULT = [(KmeansGcContentAnnotator.HIGH_GC_CONTENT_LABEL, 0, 15),
                           (KmeansGcContentAnnotator.LOW_GC_CONTENT_LABEL, 15, 32),
                           (KmeansGcContentAnnotator.HIGH_GC_CONTENT_LABEL, 32, 34),
                           (KmeansGcContentAnnotator.LOW_GC_CONTENT_LABEL, 34, 46),
                           (KmeansGcContentAnnotator.HIGH_GC_CONTENT_LABEL, 46, 67)]

        annotator = KmeansGcContentAnnotator(kmer_length=3, minimum_feature_size=None)

        annotation_labels = annotator.annotate_sequence(INPUT_SEQUENCE)
        annotation_labels = [self.intervalToTupple(interval) for interval in annotation_labels]

        self.assertListEqual(annotation_labels, EXPECTED_RESULT)

    def testEasyExampleWithSmoothing(self):
        """
        Test the states with smoothing
        """

        INPUT_SEQUENCE = "GCGCCCCGCAGCGCGATATATATATATAATATGCATATATATATATGCGCGCGCGCGGCGCGCGCGC"

        EXPECTED_RESULT = [(KmeansGcContentAnnotator.HIGH_GC_CONTENT_LABEL, 0, 15),
                           (KmeansGcContentAnnotator.LOW_GC_CONTENT_LABEL, 15,  46),
                           (KmeansGcContentAnnotator.HIGH_GC_CONTENT_LABEL, 46, 67)]

        annotator = KmeansGcContentAnnotator(kmer_length=3, minimum_feature_size=5)

        annotation_labels = annotator.annotate_sequence(INPUT_SEQUENCE)
        annotation_labels = [self.intervalToTupple(interval) for interval in annotation_labels]

        self.assertListEqual(annotation_labels, EXPECTED_RESULT)

    def testEasyExampleWithSmoothingGCOnly(self):
        """
        Test where there is only one state throughout the whole sequence
        """

        INPUT_SEQUENCE = "GCGCCCCGCGCGCGGCGCGGCGCGCGCGC"

        EXPECTED_RESULT = [StateInterval(KmeansGcContentAnnotator.UNKNOWN_LABEL, 0, 29)]

        annotator = KmeansGcContentAnnotator(kmer_length=3, minimum_feature_size=5)

        annotation_labels = annotator.annotate_sequence(INPUT_SEQUENCE)

        self.assertListEqual(annotation_labels, EXPECTED_RESULT)

    def testScoringNotEnoughRegions(self):
        """
        Test scoring with an obvious region of GC bias
        """

        INPUT_SEQUENCE = "GCGCCCCGCAGCTAGGCGCGGCGCGATATATATATATAATATGCATATATATATATGCAGCGCGCGCGCGCGGCGCGCGCGC"

        annotator = KmeansGcContentAnnotator(kmer_length=10, minimum_feature_size=5)

        annotation_labels = annotator.annotate_sequence(INPUT_SEQUENCE)

        for interval in annotation_labels:
            msg = "Interval {}: should not have enough data to give a pval".format(interval)
            self.assertIsNone(interval.score, msg)

    def testScoringDiscenibleHighGCRegion(self):
        """
        Test scoring with a highly Discernible GC region
        :return:
        """
        INPUT_SEQUENCE = "GTGCACCACAGCTAGGCACGGCTATGCGCATGCATGCGAATGCGGCATCTCGAGGGCCATGCATATATATATATATATATATATATATATATATATATATATATGCAGCGCGCGCGCGCGGCGCGCGCGC"
        ALPHA = 0.05

        annotator = KmeansGcContentAnnotator(kmer_length=10, minimum_feature_size=10)

        annotation_labels = annotator.annotate_sequence(INPUT_SEQUENCE)
        condition_passed = False
        msg = "could not find the region with high GC content in these elements: {}".format(annotation_labels)
        for interval in annotation_labels:
            if interval.start == 104:
                condition_passed = interval.score < ALPHA

        self.assertTrue(condition_passed, msg)

if __name__ == "__main__":
    unittest.main()
