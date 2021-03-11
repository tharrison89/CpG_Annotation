"""
An annotator that annotates based on relative GC content
"""

from annotator.annotator_base import AnnotatorBase, StateInterval
import annotator.annotation_utils as utils
from statsmodels.stats.proportion import proportions_ztest

class GCScoredInterval(StateInterval):
    """
    An interval specific to GC content applications
    """

    def __init__(self, label, start, stop):
        """
        ctor

        :param str label: State label
        :param int start: position start (expects a half open interval)
        :param int stop: position stop (expects a half open interval)
        """

        super(GCScoredInterval, self).__init__(label, start, stop)
        self._gc_content = None
        self._at_bases = None
        self._gc_bases = None

    def set_sequence_composition_stats(self, sequence):
        """
        Given the contig compute GC stats over the interval

        :param str sequence: The contig
        """
        subsequence = sequence[self.start: self.stop]
        self._gc_content = utils.compute_gc_content(subsequence, len(subsequence))[0][0]
        self._gc_bases = round(len(subsequence) * self._gc_content)
        self._at_bases = round(len(subsequence) - self._gc_bases)

    def __str__(self):
        gc_content = "UNKNOWN_GC_CONTENT"
        if self._gc_content is not None:
            gc_content = self._gc_content

        if self.score is None:
            return "{}\tNot Enough Data\t{}".format(super(GCScoredInterval, self).__str__(), gc_content)
        else:
            return "{}\t{}\t{}".format(super(GCScoredInterval, self).__str__(), self.score, gc_content)

    def __repr__(self):
        return self.__str__()

    @property
    def gc_content(self):
        """
        GC content property
        """
        return self._gc_content

    @property
    def gc_count(self):
        """
        GC count property
        """
        return self._gc_bases

    @property
    def at_count(self):
        """
        AT counts
        """
        return self._at_bases


class KmeansGcContentAnnotator(AnnotatorBase):
    """
    An annotator that will anotate in two states based on relative GC content
    """
    DEFAULT_MIN_FEATURE_SIZE = 10
    DEFAULT_KMER_LENGTH = 50
    _NUMBER_OF_STATES = 2
    MIN_REQUIRED_BASELINE_BASES = 30
    DEFAULT_ALPHA = 0.01

    HIGH_GC_CONTENT_LABEL = "High GC Content"
    LOW_GC_CONTENT_LABEL = "Low GC Content"
    UNKNOWN_LABEL = "UNKNOWN"

    def __init__(self, minimum_feature_size=DEFAULT_MIN_FEATURE_SIZE, kmer_length=DEFAULT_KMER_LENGTH):
        """
        ctor

        :param int minimum_feature_size: Minimum feature size for state smoothing
        :param int kmer_length: Kmer length for sliding window
        """
        self._kmer_length = kmer_length
        self._label_mapping = None
        super(KmeansGcContentAnnotator, self).__init__(number_of_states=self._NUMBER_OF_STATES,
                                                       minimum_feature_size=minimum_feature_size)
        self._interval_type = GCScoredInterval

    def _compute_state_sequence(self, input_sequence):
        """
        Annotate the states of sequence

        :param str input_sequence: Input sequence
        :return list(int): List of sequence states for each position in the sequence
        """

        self._label_mapping = {}
        gc_content_sequence = utils.compute_gc_content(input_sequence, kmer_size=self._kmer_length)
        states = utils.kmeans(gc_content_sequence, number_of_clusters=self._NUMBER_OF_STATES)

        first_state = states[0]
        last_state = states[-1]

        # we want GC content to be in the middle of the window
        for i in xrange(self._kmer_length / 2):
            states.insert(0, first_state)



        # since this is a sliding window approach there will be a coupled unlabeled position just
        # continue with the final state
        while len(states) < len(input_sequence):
            states.append(last_state)

        candidates = {}
        for i in xrange(len(states)):
            if len(candidates) == 2:
                break
            if not states[i] in candidates:
                candidates[states[i]] = gc_content_sequence[i][0]


        # find the kmers with the highest GC content and cache it for later use when the intervals are formed
        # This will be much quicker then labeling each state
        keys = candidates.keys()
        if len(candidates) == 1:
            self._label_mapping[keys[0]] = self.UNKNOWN_LABEL

        elif len(candidates) == 2:
            if candidates[keys[0]] > candidates[keys[1]]:
                self._label_mapping[keys[0]] = self.HIGH_GC_CONTENT_LABEL
                self._label_mapping[keys[1]] = self.LOW_GC_CONTENT_LABEL
            else:
                self._label_mapping[keys[0]] = self.LOW_GC_CONTENT_LABEL
                self._label_mapping[keys[1]] = self.HIGH_GC_CONTENT_LABEL

        return states

    def _compute_annotation_labels(self, state_sequence):
        """
        This will do the smoothing found in the base class and will the relabel clusters based on relative GC content
        :param list(int) state_sequence: The sequence of labeled states
        :return list(StateInterval): State annotations
        """
        intervals = super(KmeansGcContentAnnotator, self)._compute_annotation_labels(state_sequence)
        for interval in intervals:
            interval.relabel(self._label_mapping[interval.label])

        self._annotation = intervals
        return self._annotation

    def _assign_score(self, annotation_labels, input_sequence):
        """
        Assigns a p-value to the interval based on how far the GC content is according to the opposite state

        :param list(StateInterval) annotation_labels: State intervals
        :param str input_sequence: The input sequence
        """
        # Computes total GC content, A-T (and N), and GC count for each region
        for annotation_interval in annotation_labels:
            annotation_interval.set_sequence_composition_stats(input_sequence)

        gc_intervals = []
        at_intervals = []

        for annotation_interval in annotation_labels:
            if annotation_interval.label == self.HIGH_GC_CONTENT_LABEL:
                gc_intervals.append(annotation_interval)

            elif annotation_interval.label == self.LOW_GC_CONTENT_LABEL:
                at_intervals.append(annotation_interval)

        gc_baseline_at_counts, gc_baseline_gc_counts = self._get_basline_counts(gc_intervals)

        at_baseline_at_counts, at_baseline_gc_counts = self._get_basline_counts(at_intervals)

        for interval in annotation_labels:
            # observed is the counts of GC's vs. ATNs we see in the interval
            observed_counts = [interval.at_count, interval.gc_count]
            interval.score = None

            # While the baseline is pulled from the opposite state
            if interval.label == self.HIGH_GC_CONTENT_LABEL:
                if gc_baseline_at_counts is None:
                    continue
                expected_counts = [gc_baseline_at_counts, gc_baseline_gc_counts]
                alternative_hypothesis = "smaller"

            elif interval.label == self.LOW_GC_CONTENT_LABEL:
                if gc_baseline_at_counts is None:
                    continue
                expected_counts = [at_baseline_at_counts, at_baseline_gc_counts]
                alternative_hypothesis = "larger"
            else:
                continue

            expected_counts[0] -= interval.at_count
            expected_counts[1] -= interval.gc_count

            if sum(expected_counts) < self.MIN_REQUIRED_BASELINE_BASES:
                continue

            observed_true_trials = observed_counts[0]
            observed_trials = sum(observed_counts)

            expected_true_trials = expected_counts[0]
            expected_trials = sum(expected_counts)

            true_trials = [observed_true_trials, expected_true_trials]

            depths = [observed_trials, expected_trials]

            interval.score = proportions_ztest(count=true_trials, nobs=depths,  alternative=alternative_hypothesis)[1]

    @classmethod
    def _get_basline_counts(cls, intervals):
        """
        Helper method to get baseline counts for the scoring test

        :param list(StateInterval) annotation_labels: State intervals with GC content info filled out
            sorted by priority to add to the baseline
        :return int/None, int/None: AT vs GC counts to go in the baseline or Nones if there were not enough
            info to generate a baseline
        """

        at_counts = 0
        gc_counts = 0

        for interval in intervals:
            at_counts += interval.at_count
            gc_counts += interval.gc_count


        return at_counts, gc_counts

    def iterate_significant_hits(self, alpha=DEFAULT_ALPHA):
        """
        Iterate significant hits given an alpha

        :param float alpha: Alpha for significance
        :return genrator(GCScoredInterval): Significant intervals
        """
        return (interval for interval in self._annotation if interval.score <= alpha)
