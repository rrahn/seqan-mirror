/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
  ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ============================================================================
  Read Analysis code for 454 reads.
  ==========================================================================*/

#ifndef READ_ANALYZER_ANALYZE_454_H_
#define READ_ANALYZER_ANALYZE_454_H_

#include "read_analyzer.h"
#include "histogram.h"

using namespace seqan;

// ============================================================================
// Enums, Classes
// ============================================================================

struct LS454_ {};
typedef Tag<LS454_> LS454;

template<>
struct ReadEvaluationResult<LS454>
{
public:
    // Histogram over the length of reads.
    Histogram<Dense> readLengths;
    // Histogram over the number of homopolymers in each read.
    Histogram<Dense> homopolymerCounts;
};


template<>
struct AlignmentEvaluationResult<LS454>
{
public:
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TFragmentStore>
void performReadEvaluation(ReadEvaluationResult<LS454> & result, TFragmentStore & fragmentStore)
{
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename Iterator<TReadSeqStore, Standard>::Type TReadSeqStoreIterator;
    typedef typename Value<TReadSeqStoreIterator>::Type TString;
    typedef typename Value<TString>::Type TAlphabet;

    // Create histograms of the read length and number of homopolymers.
    for (TReadSeqStoreIterator it = begin(fragmentStore.readSeqStore); it != end(fragmentStore.readSeqStore); ++it) {
        tally(result.readLengths, length(*it));

        unsigned homopolymerCount = 1;
        TAlphabet c = (*it)[0];
        for (unsigned i = 0; i < length(*it); ++i) {
            if ((*it)[i] != c) {
                c = (*it)[i];
                homopolymerCount += 1;
            }
        }
        tally(result.homopolymerCounts, homopolymerCount);
    }
}


void printReadEvaluationResults(ReadEvaluationResult<LS454> const & result)
{
    std::cout << "#--file:read-lengths.dat" << std::endl;
    std::cout << "#Read Length  Occurences" << std::endl;
    for (unsigned i = 0; i < length(result.readLengths.counters); ++i)
        if (result.readLengths.counters[i] > 0u)
            printf("%12u  %10.f\n", i, result.readLengths.counters[i]);

    std::cout << "#--file:homopolymer-counts.dat" << std::endl;
    std::cout << "#Homopolymer Count  Occurences" << std::endl;
    for (unsigned i = 0; i < length(result.homopolymerCounts.counters); ++i)
        if (result.homopolymerCounts.counters[i] > 0u)
            printf("%12u  %10.f\n", i, result.homopolymerCounts.counters[i]);
}

#endif  // #ifndef READ_ANALYZER_ANALYZE_454_H_
