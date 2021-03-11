"""
Base class for annotating a sequence with a state
"""

import abc
import six

class StateInterval(object):
    """
    Object representing a spanning region where everything is one state
    """

    def __init__(self, label, start, stop):
        """
        ctor

        :param str label: State label
        :param int start: position start (expects a half open interval)
        :param int stop: position stop (expects a half open interval)
        """
        self._label = label
        self._start = start
        self._stop = stop
        self.score = None

    def __iter__(self):
        """
        Magic method for iteration
        :yields int: positions
        """
        return xrange(self.start, self.stop)

    def __str__(self):
        """
        Magic method to represent as a string

        :return str: self in stringified form
        """
        return "{}: {}-{}".format(self._label, self._start, self._stop)

    def __repr__(self):
        """
        Magic method to represent as a string

        :return str: self in stringified form
        """
        return self.__str__()

    def __eq__(self, other):
        """
        Magic method for equality

        :param StateInterval other: The other being compared
        :return bool: Are these equivolent?
        """
        return str(self._label) == str(other.label) and self.start == other.start and self.stop == other.stop

    def relabel(self, new_label):
        """
        Relabel the state

        :param str new_label: Desired new label
        """
        self._label = new_label

    @property
    def label(self):
        """
        The label property
        """
        return self._label

    @property
    def start(self):
        """
        The start property
        """
        return self._start

    @property
    def stop(self):
        """
        The stop property
        """
        return self._stop

    def __len__(self):
        """
        Magic method for returning a length\

        :return int: Length of interval
        """
        return self._stop - self._start

@six.add_metaclass(abc.ABCMeta)
class AnnotatorBase(object):
    """
    Base class for annotating a sequence with states
    """

    def __init__(self, number_of_states, minimum_feature_size=None):
        """
        ctor

        :param int number_of_states: Desired number of states
        :param int minimum_feature_size: Smallest size an interval can be for smoothing
        """
        self._number_of_states = number_of_states
        self._minimum_feature_size = minimum_feature_size
        self._annotation = None
        self._interval_type = StateInterval


    def annotate_sequence(self, input_sequence):
        """
        anotate a sequence

        :param str input_sequence: Input sequence for annotating
        :return list(StateInterval): State annotations
        """
        state_sequence = self._compute_state_sequence(input_sequence)
        annotation_labels = self._compute_annotation_labels(state_sequence)
        self._assign_score(annotation_labels, input_sequence)
        return annotation_labels

    @abc.abstractmethod
    def _compute_state_sequence(self, input_sequence):
        """
        Abstract method for annotating the states of sequence

        :param str input_sequence: Input sequence
        :return list(object): List of sequence states for each position in the sequence
        """
        pass

    @abc.abstractmethod
    def _assign_score(self, annotation_labels, input_sequence):
        """
        Abstract method for scoring a set of annotation labels

        :param iterable(StateInterval) annotation_labels: Intervals to score
        :param str input_sequence: Input sequence that is currently being scored
        """
        pass

    def _compute_annotation_labels(self, state_sequence):
        """
        Cluster the sequence together and provide smoothing if desired

        :param list(object) state_sequence: State sequences
        :return list(StateInterval): State annotations
        """

        intervals = []

        start = None
        prev_input = None

        for i, state in enumerate(state_sequence):
            if state != prev_input:
                if not prev_input is None:
                    stop = i
                    interval = self._interval_type(prev_input, start, stop)
                    intervals.append(interval)

                start = i
            prev_input = state

        stop = len(state_sequence)
        interval = self._interval_type(prev_input, start, stop)
        intervals.append(interval)

        if self._minimum_feature_size is not None:
            self._smooth_intervals(intervals)

        self._annotation = intervals
        return self._annotation

    def _smooth_intervals(self, intervals):
        """
        Smooth out state sequence intervals

        :param list(StateInterval) intervals: Intervals to be smoothed note this is modified directly
        """
        i = 0
        while i + 1 < len(intervals):
            interval_one = intervals[i]
            interval_one_length = len(interval_one)

            interval_two = intervals[i + 1]
            interval_two_length = len(interval_two)

            # Case adjacent intervals are sufficient size and are a different state
            if interval_one_length > self._minimum_feature_size and interval_two_length > self._minimum_feature_size \
                and interval_one.label != interval_two.label:
                i += 1
                continue
            else:
                # These two need to be combined find the smaller interval
                label_to_replace = interval_one.label \
                    if interval_two_length < interval_one_length else interval_two.label

                # combine them to create a new one than insert in place
                new_interval = self._interval_type(label=label_to_replace, start=interval_one.start,
                                                   stop=interval_two.stop)
                intervals.pop(i + 1)
                intervals.pop(i)
                intervals.insert(i, new_interval)

                # we need to go back one now in case the previous interval is the same state
                i = max(0, i - 1)
