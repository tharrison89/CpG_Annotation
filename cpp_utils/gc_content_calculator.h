/*
 * gccontentenumerator.h
 *
 *  Created on: Feb 14, 2021
 *      Author: thomas
 */

#ifndef GC_CONTENT_CALCULATOR_H_
#define GC_CONTENT_CALCULATOR_H_

#include <unordered_set>

#include "../cpp_utils/kmer_calculator.h"

/*
 * Class to compute a sliding window GC content of a string Not is of type Kmer calculator
 */
class GcContentCalculator: public KmerCalculator {
public:
	/*
	 * ctor
	 *
	 * @param kmer_size desired kmer size for sliding window
	 */
	GcContentCalculator(const int kmer_size);

    /*
     * dtor
     */
	virtual ~GcContentCalculator();

private:

	// expected GC symbols
	const std::unordered_set<char> GC_SYMBOLS{'G', 'C', 'g', 'c'};

	/*
	 * Compute metrics in this case it is one-dimensional GC content per sliding window overides base class
	 *
	 * @param kmer The input sequence
	 * @return The GC content of the kmer
	 */
	virtual std::vector<double> ComputeMetrics(const std::string kmer);

	/*
	 * Overides to speed up the computation on succesive GC calculations
	 *
	 * @param kmer Input sequence
	 * @param prev the previous vector of metrics
	 * @param prev_nuc The previous first character of the prior kmer
	 * @return The GC content of the kmer
	 */
	virtual std::vector<double> ComputeMetrics(const std::string kmer, const std::vector<double> prev, const char prev_nuc);
};

#endif /* GC_CONTENT_CALCULATOR_H_ */
