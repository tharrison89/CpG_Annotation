/*
 * General CPP utils to expose to the python class
 */

#include <iostream>
#include <string>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../cpp_utils/gc_content_calculator.h"

/*
 * Wrapper around the GC content calculator
 *
 * @param input The input sequence
 * @param kmer_size The size of the sliding window
 * @return positional GC content annotations
 */
std::vector<std::vector <double>> compute_gc_content(const std::string input, const int kmer_size){
	GcContentCalculator calculator(kmer_size);
	return calculator.KmerizeAndComputMetrics(input);
}

PYBIND11_MODULE(cpp_utils, m) {
    m.doc() = "CPP utils for the annotator workflow";

    m.def("compute_gc_content", &compute_gc_content, "A function which computes GC content with a sliding window");
}
