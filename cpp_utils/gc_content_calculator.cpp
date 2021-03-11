/*
 * gccontentenumerator.cpp
 *
 *  Created on: Feb 14, 2021
 *      Author: thomas
 */

#include "../cpp_utils/gc_content_calculator.h"

#include <math.h>

GcContentCalculator::GcContentCalculator(const int kmer_size): KmerCalculator{kmer_size} {
}

GcContentCalculator::~GcContentCalculator() {
}

std::vector<double> GcContentCalculator::ComputeMetrics(const std::string kmer){
	double gc_count = 0;
	std::vector<double> results;

	for (auto it = kmer.cbegin() ; it != kmer.cend(); ++it) {
	        if(this->GC_SYMBOLS.find(*it) != this->GC_SYMBOLS.end()){
	        	gc_count++;
	        }
	}
	 results.push_back(gc_count / this->kmer_size);
	 return results;
}


std::vector<double> GcContentCalculator::ComputeMetrics(const std::string kmer, const std::vector<double> prev, const char prev_nuc){
    // back calculate the number of G's and Cs from the previous calculation
	int gc_count = round(kmer.length() * prev.at(0));

	// factor out the nucleotide that will no longer be in the kmer
	if (this->GC_SYMBOLS.find(prev_nuc) != this->GC_SYMBOLS.end()){
		gc_count--;
	}

	// account for the new character
	if(this->GC_SYMBOLS.find(kmer.at(kmer.length() - 1)) != this->GC_SYMBOLS.end())
	{
		gc_count++;
	}
	std::vector<double> result;
	result.push_back(double(gc_count)/ kmer.length());

	return result;
}


