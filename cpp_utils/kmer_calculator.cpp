/*
 * kmerenumerator.cpp
 *
 *  Created on: Feb 14, 2021
 *      Author: thomas
 */

#include "../cpp_utils/kmer_calculator.h"

KmerCalculator::KmerCalculator(const int kmer_size)  : kmer_size{kmer_size} {}

KmerCalculator::~KmerCalculator() {}

std::vector<std::vector<double>> KmerCalculator::KmerizeAndComputMetrics(const std::string contig){

	std::vector<std::vector<double>>  results;
	std::vector<double> prev;
	char prev_nuc;
	for(int i=0; i + this->kmer_size <= contig.length(); i++){
		std::vector<double> curr;
		std::string kmer = contig.substr(i, this->kmer_size);
		if (prev.size() == 0){
			curr = this->ComputeMetrics(kmer);
		}

		else{
			// a chance for the child to optimize
			curr = this->ComputeMetrics(kmer, prev, prev_nuc);
		}

		results.push_back(curr);
		prev = curr;
		prev_nuc = kmer.at(0);

	}
	return results;
}


std::vector<double> KmerCalculator::ComputeMetrics(const std::string kmer, const std::vector<double> prev, const char prev_nuc){
    // This is default behavior if there is no wayto use previos calculations to inform speed up this method
	return this->ComputeMetrics(kmer);
}
