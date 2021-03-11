/*
 * kmerenumerator.h
 */

#ifndef KMER_CALCULATOR_H_
#define KMER_CALCULATOR_H_



#include <iostream>
#include <vector>
#include <string>

/*
 * Base class for performing computation on a sliding window
 */
class KmerCalculator {
public:

	// the size of kmer
	int kmer_size;

	/*
	 * ctor
	 *
	 * @param kmer_size The desired size of the sliding window
	 */
	KmerCalculator(const int kmer_size);

	/*
	 * dtor
	 */
	virtual ~KmerCalculator();

	/*
	 * Kmerize a sequence and compute any form of metrics on it
	 *
	 * @param contig The sequence to compute metric on
	 * @return a vector of vectors where the inner vector are metrics corresponding to a particular kmer
	 */
	std::vector<std::vector<double>> KmerizeAndComputMetrics(const std::string contig);

private:
	/*
	 * Compute metrics on a kmer
	 *
	 * @param kmer Input sequence
	 * @return A vector of metrics
	 */
	virtual std::vector<double> ComputeMetrics(const std::string kmer) = 0;

	/*
	 * Compute metrics on a kmer with prior knowledge
	 *     (a chance to optimize overide if this can speed up the implementation
	 *     otherwise it will default to calling the method with just a string)
	 *
	 * @param kmer Input sequence
	 * @param prev the previous vector of metrics
	 * @param prev_nuc The previous first character of the prior kmer
	 */
	virtual std::vector<double> ComputeMetrics(const std::string kmer, const std::vector<double> prev, const char prev_nuc);

};

#endif /* KMER_CALCULATOR_H_ */
