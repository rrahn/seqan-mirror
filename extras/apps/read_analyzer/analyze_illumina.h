/*====================================mismatchCountsPerMismatchPerPosition======================================
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
  Read Analysis code for Illumina reads.
  ==========================================================================*/

#ifndef READ_ANALYZER_ANALYZE_ILLUMINA_H_
#define READ_ANALYZER_ANALYZE_ILLUMINA_H_

#include "read_analyzer.h"

using namespace seqan;

// ============================================================================
// Enums, Classes
// ============================================================================

struct Illumina_ {};
typedef Tag<Illumina_> Illumina;


template <>
struct ReadEvaluationResult<Illumina>
{
public:
    String<size_t> baseCountOverall;  // arr[base]
    String<String<size_t> > baseCountPerPosition;  // arr[base][pos]
    String<String<String<size_t> > > qualityCountsPerPositionAndBase;  // arr[base][quality][pos]

    ReadEvaluationResult() {}
};


template <>
struct AlignmentEvaluationResult<Illumina>
{
public:
    // Note: "Mismatch" also stores matches.

    String<double> insertCountsPerBase;  // arr[base]
    String<double> deleteCountsPerBase;  // arr[base]
    String<double> mismatchCountsPerMismatch;  // arr[source base * 5 + target base]

    String<String<double> > qualityCountsForInsertPerBase;  // arr[base][quality]
    String<String<double> > qualityCountsForMismatchPerBase;  // arr[src*5+tgt][quality]

    String<String<double> > insertCountsPerBasePerPosition;  // arr[base][pos]
    String<String<double> > deleteCountsPerBasePerPosition;  // arr[base][pos]
    String<String<double> > mismatchCountsPerMismatchPerPosition;  // arr[src*5+tgt][pos]

    String<String<String<double> > > qualityCountsForInsertPerBasePerPosition;  // arr[base][quality][pos]
    String<String<String<double> > > qualityCountsForMismatchPerMismatchPerPosition;  // arr[src*5+tgt][quality][pos]
    String<String<String<double> > > qualityCountsBeforeDeletesPerBasePerPosition;  // arr[base][quality][pos]
    String<String<String<double> > > qualityCountsAfterDeletesPerBasePerPosition;  // arr[base][quality][pos]

    String<double> readsOnContig;
    String<size_t> alignmentsOnContig;

    String<double> readsWithErrors;
    String<size_t> alignmentsWithErrors;

    AlignmentEvaluationResult() {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/*
.Function.setReadLength
..summary:Set the maximal read lenth for an ReadEvaluationResult<Illumina> object.  This adjusts the internal buffer sizes and has to be called before using the object.
*/
void setReadLength(ReadEvaluationResult<Illumina> & result, unsigned readLength)
{
    clear(result.baseCountOverall);
    resize(result.baseCountOverall, 5, 0);

    clear(result.baseCountPerPosition);
    resize(result.baseCountPerPosition, 5);
    for (unsigned i = 0; i < 5; ++i)
        resize(result.baseCountPerPosition[i], readLength, 0);

    clear(result.qualityCountsPerPositionAndBase);
    resize(result.qualityCountsPerPositionAndBase, 5);
    for (unsigned i = 0; i < 5; ++i) {
        resize(result.qualityCountsPerPositionAndBase[i], 63);
        for (unsigned j = 0; j < 63; ++j) {
            resize(result.qualityCountsPerPositionAndBase[i][j], readLength, 0);
        }
    }
}

/*
.Function.countBaseWithQualityAtPosition
..summary:Tally the given (quality, base, position) combination in the EvaluationResult object.
*/
inline void countBaseWithQualityAtPosition(ReadEvaluationResult<Illumina> & result, Dna5Q base, size_t position)
{
    // Count overall bases.
    result.baseCountOverall[ordValue(base)] += 1;

    // Count overall bases per position.
    result.baseCountPerPosition[ordValue(base)][position] += 1;

    // Count overall quality per base per position.
    int b = ordValue(base);
    int quality = getQualityValue(base);
    result.qualityCountsPerPositionAndBase[b][quality][position] += 1;
}

/*
.Function.performReadEvaluation
..summary:Perform an evaluation of the read base counts and qualities.
*/
template <typename TFragmentStore>
void performReadEvaluation(ReadEvaluationResult<Illumina> & result, TFragmentStore & fragmentStore)
{
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename Iterator<TReadSeqStore, Standard>::Type TReadSeqStoreIterator;

    for (TReadSeqStoreIterator it = begin(fragmentStore.readSeqStore); it != end(fragmentStore.readSeqStore); ++it) {
        for (unsigned i = 0; i < length(*it); ++i)
            countBaseWithQualityAtPosition(result, (*it)[i], i);
    }
}

/*
.Function.printReadEvaluationResult<Illumina>s
..summary:Print statistical metrics about the read evaluation results.
*/
void printReadEvaluationResults(ReadEvaluationResult<Illumina> const & result)
{
    std::cout << "#--file:base-frequencies.dat" << std::endl;
    std::cout << "#Overall Base Frequencies" << std::endl;
    size_t sum = 0;
    for (unsigned i = 0; i < 5; ++i)
        sum += result.baseCountOverall[i];
    std::cout << "#base  ratio [%]      count" << std::endl;
    for (unsigned i = 0; i < 5; ++i) {
        std::cout << "    " << Dna5(i) << " ";
        printf("%10.2f %9lu\n", 100.0 * result.baseCountOverall[i] / sum, static_cast<long unsigned>(result.baseCountOverall[i]));
    }

    std::cout << std::endl << std::endl << "#--file:base-frequencies-position.dat" << std::endl;
    std::cout << "#Base Frequencies [%] Per Position" << std::endl;
    std::cout << "#position      A      C      G      T      N" << std::endl;
    for (unsigned i = 0; i < length(result.baseCountPerPosition[0]); ++i) {  // position
        size_t sum = 0;
        for (unsigned j = 0; j < 5; ++j)  // base
            sum += result.baseCountPerPosition[j][i];
        printf("     %4u %6.2f %6.2f %6.2f %6.2f %6.2f\n",
               i,
               100.0 * result.baseCountPerPosition[0][i] / sum,
               100.0 * result.baseCountPerPosition[1][i] / sum,
               100.0 * result.baseCountPerPosition[2][i] / sum,
               100.0 * result.baseCountPerPosition[3][i] / sum,
               100.0 * result.baseCountPerPosition[4][i] / sum);
    }

    std::cout << std::endl << std::endl << "#--file:qualities-base-position.dat" << std::endl;
    std::cout << "#Mean Quality/Std Dev Per Base Per Position" << std::endl;
    std::cout << "#position      A   sd A      C   sd C      G   sd G      T   sd T      N   sd N      *   sd *" << std::endl;
    for (unsigned i = 0; i < length(result.baseCountPerPosition[0]); ++i) {  // position
        printf("     %4u", i);
        for (unsigned j = 0; j < 5; ++j) {  // base
            // Compute mean.
            size_t sum = 0;
            size_t count = 0;
            for (unsigned k = 0; k < 63; ++k) {  // qualities
                count += result.qualityCountsPerPositionAndBase[j][k][i];
                sum += k * result.qualityCountsPerPositionAndBase[j][k][i];
            }
            double mean = 1.0 * sum / count;
            if (count == 0) {
                printf(" %6s %6s", "-", "-");
            } else {
                // Compute standard deviation.
                double devSum = 0;
                for (unsigned k = 0; k < 63; ++k) {  // qualities
                    double x = k - mean;
                    devSum += result.qualityCountsPerPositionAndBase[j][k][i] * x * x;
                }
                double stdDev = sqrt(devSum / count);
                printf(" %6.2f %6.2f", mean, stdDev);
            }
        }
        {
            // Compute mean.
            size_t sum = 0;
            size_t count = 0;
            for (unsigned j = 0; j < 5; ++j) {  // base
                for (unsigned k = 0; k < 63; ++k) {  // qualities
                    count += result.qualityCountsPerPositionAndBase[j][k][i];
                    sum += k * result.qualityCountsPerPositionAndBase[j][k][i];
                }
            }
            double mean = 1.0 * sum / count;
            if (count == 0) {
                printf(" %6s %6s", "-", "-");
            } else {
                // Compute standard deviation.
                double devSum = 0;
                for (unsigned j = 0; j < 5; ++j) {  // base
                    for (unsigned k = 0; k < 63; ++k) {  // qualities
                        double x = k - mean;
                        devSum += result.qualityCountsPerPositionAndBase[j][k][i] * x * x;
                    }
                }
                double stdDev = sqrt(devSum / count);
                printf(" %6.2f %6.2f", mean, stdDev);
            }
        }
        printf("\n");
    }
}

/*
.Function.setReadLength
..summary:Set the maximal read lenth for an AlignmentEvaluationResult<Illumina> object.  This adjusts the internal buffer sizes and has to be called before using the object.
*/
void setReadLength(AlignmentEvaluationResult<Illumina> & result, unsigned readLength, unsigned contigCount)
{
    clear(result.insertCountsPerBase);
    resize(result.insertCountsPerBase, 5, 0);
    clear(result.deleteCountsPerBase);
    resize(result.deleteCountsPerBase, 5, 0);
    clear(result.mismatchCountsPerMismatch);
    resize(result.mismatchCountsPerMismatch, 25, 0);

    clear(result.qualityCountsForInsertPerBase);
    resize(result.qualityCountsForInsertPerBase, 5);
    for (int i = 0; i < 5; ++i)
        resize(result.qualityCountsForInsertPerBase[i], 63, 0);
    clear(result.qualityCountsForMismatchPerBase);
    resize(result.qualityCountsForMismatchPerBase, 25);
    for (int i = 0; i < 25; ++i)
        resize(result.qualityCountsForMismatchPerBase[i], 63, 0);
    
    clear(result.insertCountsPerBasePerPosition);
    resize(result.insertCountsPerBasePerPosition, 5);
    for (int i = 0; i < 5; ++i)
        resize(result.insertCountsPerBasePerPosition[i], readLength, 0);
    clear(result.deleteCountsPerBasePerPosition);
    resize(result.deleteCountsPerBasePerPosition, 5);
    for (int i = 0; i < 5; ++i)
        resize(result.deleteCountsPerBasePerPosition[i], readLength, 0);
    clear(result.mismatchCountsPerMismatchPerPosition);
    resize(result.mismatchCountsPerMismatchPerPosition, 25);
    for (int i = 0; i < 25; ++i)
        resize(result.mismatchCountsPerMismatchPerPosition[i], readLength, 0);

    clear(result.qualityCountsBeforeDeletesPerBasePerPosition);
    resize(result.qualityCountsBeforeDeletesPerBasePerPosition, 5);
    for (int i = 0; i < 5; ++i) {
        resize(result.qualityCountsBeforeDeletesPerBasePerPosition[i], 63);
        for (int j = 0; j < 63; ++j)
            resize(result.qualityCountsBeforeDeletesPerBasePerPosition[i][j], readLength, 0);
    }
    clear(result.qualityCountsAfterDeletesPerBasePerPosition);
    resize(result.qualityCountsAfterDeletesPerBasePerPosition, 5);
    for (int i = 0; i < 5; ++i) {
        resize(result.qualityCountsAfterDeletesPerBasePerPosition[i], 63);
        for (int j = 0; j < 63; ++j)
            resize(result.qualityCountsAfterDeletesPerBasePerPosition[i][j], readLength, 0);
    }
    clear(result.qualityCountsForInsertPerBasePerPosition);
    resize(result.qualityCountsForInsertPerBasePerPosition, 5);
    for (int i = 0; i < 5; ++i) {
        resize(result.qualityCountsForInsertPerBasePerPosition[i], 63);
        for (int j = 0; j < 63; ++j)
            resize(result.qualityCountsForInsertPerBasePerPosition[i][j], readLength, 0);
    }
    clear(result.qualityCountsForInsertPerBasePerPosition);
    resize(result.qualityCountsForInsertPerBasePerPosition, 5);
    for (int i = 0; i < 5; ++i) {
        resize(result.qualityCountsForInsertPerBasePerPosition[i], 63);
        for (int j = 0; j < 63; ++j)
            resize(result.qualityCountsForInsertPerBasePerPosition[i][j], readLength, 0);
    }

    clear(result.qualityCountsForMismatchPerMismatchPerPosition);
    resize(result.qualityCountsForMismatchPerMismatchPerPosition, 25);
    for (int i = 0; i < 25; ++i) {
        resize(result.qualityCountsForMismatchPerMismatchPerPosition[i], 63);
        for (int j = 0; j < 63; ++j)
            resize(result.qualityCountsForMismatchPerMismatchPerPosition[i][j], readLength, 0);
    }

    clear(result.readsOnContig);
    resize(result.readsOnContig, contigCount, 0);
    clear(result.alignmentsOnContig);
    resize(result.alignmentsOnContig, contigCount, 0);

    clear(result.readsWithErrors);
    resize(result.readsWithErrors, 2 * readLength, 0);
    clear(result.alignmentsWithErrors);
    resize(result.alignmentsWithErrors, 2 * readLength, 0);
}

inline void
countInsertAtPositionWithBase(AlignmentEvaluationResult<Illumina> & result,
                              size_t position,
                              Dna5Q readBase,
                              unsigned readAlignmentCount)
{
    int b = ordValue(readBase);
    int q = getQualityValue(readBase);

    double x = 1.0 / readAlignmentCount;

    result.insertCountsPerBase[b] += x;
    result.qualityCountsForInsertPerBase[b][q] += x;
    SEQAN_ASSERT_LT(position, length(result.insertCountsPerBasePerPosition[b]));
    result.insertCountsPerBasePerPosition[b][position] += x;
    result.qualityCountsForInsertPerBasePerPosition[b][q][position] += x;
}

inline void
countDeleteAtPositionWithBase(AlignmentEvaluationResult<Illumina> & result,
                              size_t position,
                              Dna5 referenceBase,
                              bool hasBefore,
                              Dna5Q beforeReadBase,
                              bool hasAfter,
                              Dna5Q afterReadBase,
                              unsigned readAlignmentCount)
{
    int qb = getQualityValue(beforeReadBase);
    int qa = getQualityValue(afterReadBase);

    int b = ordValue(referenceBase);

    double x = 1.0 / readAlignmentCount;

    result.deleteCountsPerBase[b] += x;
    SEQAN_ASSERT_LT(position, length(result.deleteCountsPerBasePerPosition[b]));
    result.deleteCountsPerBasePerPosition[b][position] += x;
    if (hasBefore)
        result.qualityCountsBeforeDeletesPerBasePerPosition[b][qb][position] += x;
    if (hasAfter)
        result.qualityCountsAfterDeletesPerBasePerPosition[b][qa][position] += x;
}

inline void
countMismatchAtPositionWithBase(AlignmentEvaluationResult<Illumina> & result,
                                size_t position,
                                Dna5 referenceBase,
                                Dna5Q readBase,
                                unsigned readAlignmentCount)
{
    int s = ordValue(referenceBase);
    int t = ordValue(readBase);
    int q = getQualityValue(readBase);

    // if (convert<Dna5>(readBase) != referenceBase)
    //     std::cout << position << " " << q << " " << convert<char>(readBase) << std::endl;

    double x = 1.0 / readAlignmentCount;

    result.mismatchCountsPerMismatch[s * 5 + t] += x;
    result.qualityCountsForMismatchPerBase[s * 5 + t][q] += x;
    SEQAN_ASSERT_LT(position, length(result.mismatchCountsPerMismatchPerPosition[s * 5 + t]));
    result.mismatchCountsPerMismatchPerPosition[s * 5 + t][position] += x;
    result.qualityCountsForMismatchPerMismatchPerPosition[s * 5 + t][q][position] += x;
}

template <typename TAlignedReadStoreElement>
struct LessThanReadId
{
    bool operator()(TAlignedReadStoreElement const & a, TAlignedReadStoreElement const & b) const
    {
        return a.readId < b.readId;
    }
};


// Collect the aligned reads ids of all reads with more than maxErrors
// in bogusReadIds.
template <typename TReadIdSet, typename TAlignedReadsIterator, typename TFragmentStore>
void computeBogusReads(
        TReadIdSet & bogusReadIds,
        size_t alignedReadId,
        TAlignedReadsIterator alignedReadsBegin,
        TAlignedReadsIterator alignedReadsEnd,
        TFragmentStore & fragmentStore,
        unsigned maxErrors)
{
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadsIter;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename TFragmentStore::TReadSeq TReadSeq;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef typename Value<TAlignedReadStore>::Type TAlignedReadStoreElement;
    typedef typename TAlignedReadStoreElement::TGapAnchors TReadGapAnchors;
    typedef typename TContigStoreElement::TContigSeq TContigSeq;
    typedef typename TContigStoreElement::TGapAnchors TContigGapAnchors;
    typedef Gaps<TReadSeq, AnchorGaps<TContigGapAnchors> > TReadGaps;
    typedef Gaps<TContigSeq, AnchorGaps<TContigGapAnchors> > TContigGaps;
    typedef typename Iterator<TContigGaps, Standard>::Type TContigGapAnchorsIterator;
    typedef typename Iterator<TReadGaps, Standard>::Type TReadGapAnchorsIterator;

    bogusReadIds.clear();

    for (TAlignedReadsIterator it = alignedReadsBegin; it != alignedReadsEnd; ++it, ++alignedReadId) {
        // Get contig and read sequences.
        TContigSeq const & contigSeq = fragmentStore.contigSeqStore[it->contigId].seq;
        TReadSeq readSeq = fragmentStore.readSeqStore[it->readId];
        // Get gaps for contig and read.
        TContigGaps contigGaps(contigSeq, fragmentStore.contigStore[it->contigId].gaps);
        TReadGaps readGaps(readSeq, fragmentStore.alignedReadStore[alignedReadId].gaps);
        // Limit contig gaps to aligned read position.
        setBeginPosition(contigGaps, _min(it->beginPos, it->endPos));
        setEndPosition(contigGaps, _max(it->beginPos, it->endPos));
        // Reverse-complement readSeq in-place.
        bool flipped = false;
        if (it->beginPos > it->endPos) {
            flipped = true;
            reverseComplement(readSeq);
        }

        TContigGapAnchorsIterator contigGapsIt = begin(contigGaps);
        unsigned errorCount = 0;
        for (TReadGapAnchorsIterator readGapsIt = begin(readGaps); readGapsIt != end(readGaps); ++contigGapsIt, ++readGapsIt) {
            if (isGap(readGapsIt) && isGap(contigGapsIt))
                continue;  // Skip paddings.
            if (isGap(readGapsIt)) {
                errorCount += 1;
            } else if (isGap(contigGapsIt)) {
                errorCount += 1;
            } else {
                errorCount += convert<Dna5>(*contigGapsIt) != convert<Dna5>(*readGapsIt);
            }
        }

        if (errorCount > maxErrors)
            bogusReadIds.insert(alignedReadId);
    }
}


template <typename TFragmentStore>
void performAlignmentEvaluation(AlignmentEvaluationResult<Illumina> & result, TFragmentStore & fragmentStore)
{
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadsIter;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename TFragmentStore::TReadSeq TReadSeq;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef typename Value<TAlignedReadStore>::Type TAlignedReadStoreElement;
    typedef typename TAlignedReadStoreElement::TGapAnchors TReadGapAnchors;
    typedef typename TContigStoreElement::TContigSeq TContigSeq;
    typedef typename TContigStoreElement::TGapAnchors TContigGapAnchors;
    typedef Gaps<TReadSeq, AnchorGaps<TContigGapAnchors> > TReadGaps;
    typedef Gaps<TContigSeq, AnchorGaps<TContigGapAnchors> > TContigGaps;
    typedef typename Iterator<TContigGaps, Standard>::Type TContigGapAnchorsIterator;
    typedef typename Iterator<TReadGaps, Standard>::Type TReadGapAnchorsIterator;

    // Sort aligned read by read so we can easily compute the number
    // of alignments of a read using binary search.
    sortAlignedReads(fragmentStore.alignedReadStore, SortReadId());
    // Initialize reference aligned read with a sentinel value.
    TAlignedReadStoreElement refAlignedRead;
    refAlignedRead.readId = ~0u;
    typedef std::pair<TAlignedReadsIter, TAlignedReadsIter> TAlignedReadsIterPair;
    TAlignedReadsIterPair alignedReadsBeginEnd;

    size_t alignedReadId = 0;
    unsigned readAlignmentCount = 0;
    // std::set<size_t> bogusAlignmentIds;
    for (TAlignedReadsIter it = begin(fragmentStore.alignedReadStore, Standard()); it != end(fragmentStore.alignedReadStore, Standard()); ++it, ++alignedReadId) {
        // Maybe need to re-count number of reads, if at next read id.
        if (it->readId != refAlignedRead.readId) {
            refAlignedRead.readId = it->readId;
            alignedReadsBeginEnd =
                    std::equal_range(begin(fragmentStore.alignedReadStore, Standard()),
                                     end(fragmentStore.alignedReadStore, Standard()),
                                     refAlignedRead,
                                     LessThanReadId<TAlignedReadStoreElement>());
            readAlignmentCount = alignedReadsBeginEnd.second - alignedReadsBeginEnd.first;

            // Ignore reads that have more than 10% errors.
            //unsigned maxErrors = static_cast<unsigned>(0.1 * length(fragmentStore.readSeqStore[refAlignedRead.readId]));
            // computeBogusReads(bogusAlignmentIds, alignedReadId, alignedReadsBeginEnd.first, alignedReadsBeginEnd.second, fragmentStore, maxErrors);
            // size_t bogusAlignmentCount = bogusAlignmentIds.size();
            // readAlignmentCount -= bogusAlignmentCount;
        }
        // Maybe ignore bogus read.
        // if (bogusAlignmentIds.find(alignedReadId) != bogusAlignmentIds.end())
        //     continue;
        // Count alignment / read on contig.
        result.readsOnContig[it->contigId] += 1.0 / readAlignmentCount;
        result.alignmentsOnContig[it->contigId] += 1;
        // Get contig and read sequences.
        TContigSeq & contigSeq = fragmentStore.contigStore[it->contigId].seq;
        TReadSeq readSeq = fragmentStore.readSeqStore[it->readId];
        // Get gaps for contig and read.
        TContigGaps contigGaps(contigSeq, fragmentStore.contigStore[it->contigId].gaps);
        TReadGaps readGaps(readSeq, fragmentStore.alignedReadStore[alignedReadId].gaps);
        // Limit contig gaps to aligned read position.
        setBeginPosition(contigGaps, _min(it->beginPos, it->endPos));
        setEndPosition(contigGaps, _max(it->beginPos, it->endPos));
        // Reverse-complement readSeq in-place.
        unsigned readLength = length(readSeq);
        bool flipped = false;  // TODO(holtgrew): Remove?
        if (it->beginPos > it->endPos) {
            flipped = true;
            reverseComplement(readSeq);
        }
        // std::cout << "Evaluating" << std::endl;
        // std::cout << "  " << fragmentStore.readNameStore[it->readId] << ", flipped = " << flipped << std::endl;

        unsigned errorCount = 0;
        TReadGapAnchorsIterator readGapsIt = begin(readGaps);
        TContigGapAnchorsIterator contigGapsIt = begin(contigGaps);
        unsigned readPos = 0;

        // Skip leading gaps in contig and reads. 
        //std::cout << "beginPosition(contigGaps) == " << beginPosition(contigGaps) << ", endPosition(contigGaps) == " << endPosition(contigGaps) << std::endl;
        //std::cout << "positionGapToSeq(beginPosition(contigGaps)) == " << positionGapToSeq(contigGaps, beginPosition(contigGaps)) << ", positionGapToSeq(endPosition(contigGaps)) == " << positionGapToSeq(contigGaps, endPosition(contigGaps)) << std::endl;
        //for (unsigned i = 0; i < beginPosition(row(align, 1)); ++i)
        //  ++readGapsIt;

        for (; readGapsIt != end(readGaps); ++contigGapsIt, ++readGapsIt) {
            if (isGap(readGapsIt) && isGap(contigGapsIt)) {
                //std::cout << "padding" << std::endl;
                continue;  // Skip paddings.
            }
            SEQAN_ASSERT_LT(readPos, readLength);
            unsigned reportedPos = flipped ? (readLength - readPos - 1) : readPos;
            if (isGap(readGapsIt)) {
                //std::cout << "read gap" << std::endl;
                // Deletion
                bool hasBefore = false;
                Dna5Q baseBefore('N');
                if (reportedPos > 0) {
                    TReadGapAnchorsIterator tmpIt = readGapsIt;
                    while (isGap(tmpIt) && tmpIt != begin(readGaps))
                        goPrevious(tmpIt);
                    if (!isGap(tmpIt)) {
                        hasBefore = true;
                        baseBefore = convert<Dna5Q>(*tmpIt);
                    }
                }
                bool hasAfter = false;
                Dna5Q baseAfter('N');
                if (reportedPos < readLength - 1) {
                    TReadGapAnchorsIterator tmpIt = readGapsIt;
                    while (isGap(tmpIt) && tmpIt != end(readGaps))
                        goNext(tmpIt);
                    if (tmpIt != end(readGaps)) {
                        hasAfter = true;
                        baseAfter = convert<Dna5Q>(*tmpIt);
                    }
                }
                countDeleteAtPositionWithBase(result, reportedPos, convert<Dna5Q>(*contigGapsIt), hasBefore, baseBefore, hasAfter, baseAfter, readAlignmentCount);
                errorCount += 1;
            } else if (isGap(contigGapsIt)) {
                //std::cout << "contig gap" << std::endl;
                // Insert
                countInsertAtPositionWithBase(result, reportedPos, convert<Dna5Q>(*readGapsIt), readAlignmentCount);
                readPos += 1;
                errorCount += 1;
            } else {
                //std::cout << "match/mismatch" << std::endl;
                // Match / Mismatch.
                // if (*readGapsIt != *contigGapsIt) {
                //     std::cout << "it->contigId == " << it->contigId << ", it->beginPos == " << it->beginPos << ", _it->endPos == " << it->endPos << std::endl;
                // }
                //std::cout << ">>> contig " << convert<Dna5>(*contigGapsIt) << " vs read " << convert<Dna5>(*readGapsIt) << std::endl;
                errorCount += convert<Dna5Q>(*contigGapsIt) != convert<Dna5>(*readGapsIt);
                countMismatchAtPositionWithBase(result, reportedPos, convert<Dna5Q>(*contigGapsIt), convert<Dna5Q>(*readGapsIt), readAlignmentCount);
                readPos += 1;
            }
        }
        //if (errorCount > 20) {
          //std::cout << "errors: " << errorCount << ", read: " << fragmentStore.readNameStore[it->readId] << ", contig " << fragmentStore.contigNameStore[it->contigId] << ", beginPos: " << it->beginPos << std::endl;
        //}
        if (errorCount >= length(result.readsWithErrors)) {
          std::cerr << "errorCount >= length(result.readsWithErrors)!!" << std::endl;
          exit(-1);
        }
        result.readsWithErrors[errorCount] += 1.0 / readAlignmentCount;
        result.alignmentsWithErrors[errorCount] += 1;
    }
}

template <typename TFragmentStore>
void printAlignmentEvaluationResults(AlignmentEvaluationResult<Illumina> const & result, TFragmentStore const & fragmentStore)  // TODO(holtgrew): Remove parameter fragmentStore.
{
    // Print alignments / reads on contigs.
    {
        std::cout << std::endl << std::endl << "#--file:counts-per-contig.dat" << std::endl;
        std::cout << "#contig                            reads reads% alignments alignments%" << std::endl;
        double totalReads = 0;
        size_t totalAlignments = 0;
        for (unsigned i = 0; i < length(result.readsOnContig); ++i) {
          totalAlignments += result.alignmentsOnContig[i];
          totalReads += result.readsOnContig[i];
        }
        for (unsigned i = 0; i < length(result.readsOnContig); ++i) {
          printf("%-30s %9.0f %6.2f  %9lu      %6.2f\n", toCString(fragmentStore.contigNameStore[i]),
                 result.readsOnContig[i], 100.0 * result.readsOnContig[i] / totalReads,
                 static_cast<long unsigned>(result.alignmentsOnContig[i]), 100.0 * result.alignmentsOnContig[i] / totalAlignments);
        }
        printf("%-30s %9.0f %6.2f  %9lu      %6.2f\n", "*", totalReads, 100.0, static_cast<long unsigned>(totalAlignments), 100.0);
    }

    // Print alignments / reads with errors.
    {
        std::cout << std::endl << std::endl << "#--file:reads-with-errors.dat" << std::endl;
        std::cout << "#Errors    reads reads% alignments alignments%" << std::endl;
        double totalReads = 0;
        size_t totalAlignments = 0;
        for (unsigned i = 0; i < length(result.alignmentsWithErrors); ++i) {
          totalAlignments += result.alignmentsWithErrors[i];
          totalReads += result.readsWithErrors[i];
        }
        for (unsigned i = 0; i < length(result.alignmentsWithErrors); ++i) {
          printf("%2d %9.0f %6.2f  %9lu      %6.2f\n", i,
                 result.readsWithErrors[i], 100.0 * result.readsWithErrors[i] / totalReads,
                 static_cast<long unsigned>(result.alignmentsWithErrors[i]), 100.0 * result.alignmentsWithErrors[i] / totalAlignments);
        }
        printf("%-8s %9.0f %6.2f  %9lu      %6.2f\n", "*", totalReads, 100.0, static_cast<long unsigned>(totalAlignments), 100.0);
    }

    // Print error counts per base.
    std::cout << std::endl << std::endl << "#--file:error-counts-base.dat" << std::endl;
    printf("#base       insert insert%%       delete delete%%        mismatch mismatch%%        match  match%%\n");
    double totalInserts = 0;
    double totalDeletes = 0;
    double totalMismatches = 0;
    double totalMatches = 0;
    for (int i = 0; i < 5; ++i) {
        std::cout << "    " << Dna5(i);
        double mismatches = 0;
        double matches = 0;
        for (int j = 0; j < 5; ++j)
            if (i != j)
                mismatches += result.mismatchCountsPerMismatch[i * 5 + j];
            else
                matches += result.mismatchCountsPerMismatch[i * 5 + j];
        totalMismatches += mismatches;
        totalMatches += matches;
        totalInserts += result.insertCountsPerBase[i];
        totalDeletes += result.deleteCountsPerBase[i];
        double total = result.insertCountsPerBase[i] + result.deleteCountsPerBase[i] + mismatches + matches;
        printf(" %12.2f %6.2f%% %12.2f %6.2f%%    %12.2f   %6.2f%% %12.2f %6.2f%%\n", result.insertCountsPerBase[i], 100.0 * result.insertCountsPerBase[i] / total, result.deleteCountsPerBase[i], 100.0 * result.deleteCountsPerBase[i] / total, mismatches, 100.0 * mismatches / total, matches, 100.0 * matches / total);
    }
    double total = totalInserts + totalDeletes + totalMismatches + totalMatches;
    printf("    * %12.2f %6.2f%% %12.2f %6.2f%%    %12.2f   %6.2f%% %12.2f %6.2f%%\n", totalInserts, 100.0 * totalInserts / total, totalDeletes, 100.0 * totalDeletes / total, totalMismatches, 100.0 * totalMismatches / total, totalMatches, 100.0 * totalMatches / total);

    // Print substitution counts.
    std::cout << std::endl << std::endl << "#--file:substitution-counts.dat" << std::endl;
    std::cout << "#          " << Dna5(0) << "     " << Dna5(0) << " %         " << Dna5(1) << "     " << Dna5(1) << " %        " << Dna5(2) << "      " << Dna5(2) << " %         " << Dna5(3) << "     " << Dna5(3) << " %         " << Dna5(4) << "     " << Dna5(4) << " %         *" << std::endl;
    String<double> sums;
    resize(sums, 6, 0);
    for (int i = 0; i < 5; ++i) {
        double sum = 0;
        std::cout << Dna5(i) << " ";
        for (int j = 0; j < 5; ++j) {
            sums[j] += result.mismatchCountsPerMismatch[i * 5 + j];
            sum += result.mismatchCountsPerMismatch[i * 5 + j];
        }
        for (int j = 0; j < 5; ++j) {
            if (sum > 0)
                printf(" %9.2f %6.2f%%", result.mismatchCountsPerMismatch[i * 5 + j], 100.0 * result.mismatchCountsPerMismatch[i * 5 + j] / sum);
            else
                printf(" %9.2f %6s ", result.mismatchCountsPerMismatch[i * 5 + j], "-");
        }
        printf(" %9.2f", sum);
        sums[5] += sum;
        std::cout << std::endl;
    }
    total = 0;
    for (int i = 0; i < 5; ++i)
        total += sums[i];
    std::cout << "* ";
    for (int i = 0; i < 5; ++i)
        printf(" %9.2f %6.2f%%", sums[i], 100.0 * sums[i] / total);
    printf(" %9.2f", sums[5]);
    std::cout << std::endl;

    // Print mean/stddev of qualities per base for inserts.
    std::cout << std::endl << std::endl << "#--file:insert-qualities.dat" << std::endl;
    printf("#base   mean     sd\n");
    double totalMeanSum = 0;
    double totalCount = 0;
    for (int i = 0; i < 5; ++i) {
        std::cout << "    " << Dna5(i);
        double meanSum = 0;
        double count = 0;
        for (size_t j = 0; j < 63; ++j) {
            meanSum += result.qualityCountsForInsertPerBase[i][j] * j;
            count += result.qualityCountsForInsertPerBase[i][j];
        }
        totalMeanSum += meanSum;
        totalCount += count;
        if (count > 0) {
            double mean = 1.0 * meanSum / count;
            double stdDevSum = 0;
            for (size_t j = 0; j < 63; ++j)
                stdDevSum += result.qualityCountsForInsertPerBase[i][j] * (mean - j) * (mean - j);
            double stdDev = ::std::sqrt(stdDevSum / count);
            printf(" %6.2f %6.2f", mean, stdDev);
        } else {
            printf(" %6s %6s", "-", "-");
        }
        std::cout << std::endl;
    }
    std::cout << "    *";
    if (totalCount > 0) {
        double totalMean = 1.0 * totalMeanSum / totalCount;
        double totalStdDevSum = 0;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 63; ++j)
                totalStdDevSum += result.qualityCountsForInsertPerBase[i][j] * (totalMean - j) * (totalMean - j);
        double totalStdDev = ::std::sqrt(totalStdDevSum / totalCount);
        printf(" %6.2f %6.2f", totalMean, totalStdDev);
    } else {
        printf(" %6s %6s", "-", "-");
    }
    std::cout << std::endl;

    // Print mean/stddev of qualities per base for mismatches.
    std::cout << std::endl << std::endl << "#--file:insert-mismatch.dat" << std::endl;
    printf("#base   mean     sd\n");
    totalMeanSum = 0;
    totalCount = 0;
    for (int i = 0; i < 5; ++i) {
        std::cout << "    " << Dna5(i);
        double meanSum = 0;
        double count = 0;
        for (size_t j = 0; j < 63; ++j) {
            meanSum += result.qualityCountsForInsertPerBase[i][j] * j;
            count += result.qualityCountsForInsertPerBase[i][j];
        }
        totalMeanSum += meanSum;
        totalCount += count;
        if (count > 0) {
            double mean = 1.0 * meanSum / count;
            double stdDevSum = 0;
            for (size_t j = 0; j < 63; ++j)
                stdDevSum += result.qualityCountsForInsertPerBase[i][j] * (mean - j) * (mean - j);
            double stdDev = ::std::sqrt(stdDevSum / count);
            printf(" %6.2f %6.2f", mean, stdDev);
        } else {
            printf(" %6s %6s", "-", "-");
        }
        std::cout << std::endl;
    }
    std::cout << "    *";
    if (totalCount > 0) {
        double totalMean = 1.0 * totalMeanSum / totalCount;
        double totalStdDevSum = 0;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 63; ++j)
                totalStdDevSum += result.qualityCountsForInsertPerBase[i][j] * (totalMean - j) * (totalMean - j);
        double totalStdDev = ::std::sqrt(totalStdDevSum / totalCount);
        printf(" %6.2f %6.2f", totalMean, totalStdDev);
    } else {
        printf(" %6s %6s", "-", "-");
    }
    std::cout << std::endl;

    // Qualities per substitution.
    std::cout << std::endl << std::endl << "#--file:qualities-mismatch-base.dat" << std::endl;
    std::cout << "#Mean Quality/Std Dev Per Substitution" << std::endl;
    std::cout << "# * means any, + means any but match" << std::endl;
    std::cout << "#       A   sd A      C   sd C      G   sd G      T   sd T      N   sd N      +   sd +      *   sd *" << std::endl;
    for (unsigned i = 0; i < 5; ++i) {  // source base
        std::cout << Dna5(i) << " ";
        for (unsigned j = 0; j < 5; ++j) {  // target base
            // Compute mean.
            double sum = 0;
            double count = 0;
            for (unsigned k = 0; k < 63; ++k) {  // qualities
                count += result.qualityCountsForMismatchPerBase[i * 5 + j][k];
                sum += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * k;
            }
            double mean = 1.0 * sum / count;
            if (count == 0) {
                printf(" %6s %6s", "-", "-");
            } else {
                // Compute standard deviation.
                double devSum = 0;
                for (unsigned k = 0; k < 63; ++k) {  // qualities
                    double x = k - mean;
                    devSum += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * x * x;
                }
                double stdDev = sqrt(devSum / count);
                printf(" %6.2f %6.2f", mean, stdDev);
            }
        }
        // Source to all.
        {
            // Compute mean.
            double sum = 0;
            double count = 0;
            double sumMismatch = 0;
            double countMismatch = 0;
            for (unsigned j = 0; j < 5; ++j) {  // base
                for (unsigned k = 0; k < 63; ++k) {  // qualities
                    if (i != j) {
                        sumMismatch += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * k;
                        countMismatch += result.qualityCountsForMismatchPerBase[i * 5 + j][k];
                    }
                    sum += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * k;
                    count += result.qualityCountsForMismatchPerBase[i * 5 + j][k];
                }
            }
            double mean = 1.0 * sum / count;
            double meanMismatch = 1.0 * sumMismatch / countMismatch;
            if (countMismatch == 0)
                printf(" %6s %6s", "-", "-");
            if (count == 0)
                printf(" %6s %6s", "-", "-");
            if (count != 0 || countMismatch != 0) {
                // Compute standard deviation.
                double devSum = 0;
                double devSumMismatch = 0;
                for (unsigned j = 0; j < 5; ++j) {  // base
                    for (unsigned k = 0; k < 63; ++k) {  // qualities
                        if (i != j) {
                            double x = k - meanMismatch;
                            devSumMismatch += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * x * x;
                        }
                        double x = k - mean;
                        devSum += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * x * x;
                    }
                }
                if (countMismatch != 0) {
                    double stdDevMismatch = sqrt(devSumMismatch / countMismatch);
                    printf(" %6.2f %6.2f", meanMismatch, stdDevMismatch);
                }
                if (count != 0) {
                    double stdDev = sqrt(devSum / count);
                    printf(" %6.2f %6.2f", mean, stdDev);
                }
            }
        }
        printf("\n");
    }
    // Mismatch to target.
    // TODO(holtgrew): Write me!
    std::cout << "+ TODO" << std::endl;
    // All to target.
    // TODO(holtgrew): Write me!
    std::cout << "* TODO" << std::endl;

    // Error probabilities per position.
    std::cout << std::endl << std::endl << "#--file:error-probabilities-position.dat" << std::endl;
    std::cout << "#position      insert  insert%    delete  delete%  mismatch mismatch%       any     any%     total" << std::endl;
    for (unsigned i = 0; i < length(result.insertCountsPerBasePerPosition[0]); ++i) {  // position
        double inserts = 0;
        double deletes = 0;
        double mismatches = 0;
        double matches = 0;
        for (unsigned b = 0; b < 5; ++b) {
            inserts += result.insertCountsPerBasePerPosition[b][i];
            deletes += result.deleteCountsPerBasePerPosition[b][i];
            for (unsigned a = 0; a < 5; ++a) {
                if (a == b)
                    matches += result.mismatchCountsPerMismatchPerPosition[a * 5 + b][i];
                else
                    mismatches += result.mismatchCountsPerMismatchPerPosition[a * 5 + b][i];
            }
        }
        // The number of aligned reads is fixed.  Note that there may
        // be multiple error events per position but at most one
        // mismatch.
        double total = inserts + deletes + mismatches + matches;
        printf("%9u   %9.2f %8.5f %9.2f %8.5f %9.2f  %8.5f %9.2f %8.5f %9.2f\n", i, inserts, 100.0 * inserts / total, deletes, 100.0 * deletes / total, mismatches, 100.0 * mismatches / total, inserts + deletes + mismatches, 100.0 * (inserts + deletes + mismatches) / total, total);
    }

    // Qualities per position and error type.
    std::cout << std::endl << std::endl << "#--file:qualities-per-error-position.dat" << std::endl;
    std::cout << "#position insert sd insert delete before sd delete before delete after delete after sd mismatch   sd mismatch   match sd match  total  sd total" << std::endl;
    for (unsigned position = 0; position < length(result.insertCountsPerBasePerPosition[0]); ++position) {  // position
        printf("%9u", position);
        // insert mean
        double insertSum = 0;
        double insertCount = 0;
        double insertMean = 0;
        for (int b = 0; b < 5; ++b) {
            for (int q = 0; q < 63; ++q) {
                insertSum += result.qualityCountsForInsertPerBasePerPosition[b][q][position] * q;
                insertCount += result.qualityCountsForInsertPerBasePerPosition[b][q][position];
            }
        }
        if (insertCount == 0) {
            printf("%7s   %7s", "-", "-");
        } else {
            insertMean = insertSum / insertCount;
            // insert sd
            double insertStdDevSum = 0;
            for (int b = 0; b < 5; ++b) {
                for (int q = 0; q < 63; ++q) {
                    double x = insertMean - q;
                    insertStdDevSum += result.qualityCountsForInsertPerBasePerPosition[b][q][position] * x * x;
                }
            }
            double insertStdDev = std::sqrt(1.0 / insertCount * insertStdDevSum);
            printf("%7.2f   %7.2f", insertMean, insertStdDev);
        }
        // delete mean
        double deleteBeforeSum = 0;
        double deleteBeforeCount = 0;
        double deleteBeforeMean = 0;
        double deleteAfterSum = 0;
        double deleteAfterCount = 0;
        double deleteAfterMean = 0;
        for (int b = 0; b < 5; ++b) {
            for (int q = 0; q < 63; ++q) {
                deleteBeforeSum += result.qualityCountsBeforeDeletesPerBasePerPosition[b][q][position] * q;
                deleteBeforeCount += result.qualityCountsBeforeDeletesPerBasePerPosition[b][q][position];
                deleteAfterSum += result.qualityCountsAfterDeletesPerBasePerPosition[b][q][position] * q;
                deleteAfterCount += result.qualityCountsAfterDeletesPerBasePerPosition[b][q][position];
            }
        }
        if (deleteBeforeCount == 0) {
            printf("       %7s          %7s", "-", "-");
        } else {
            deleteBeforeMean = deleteBeforeSum / deleteBeforeCount;
            // insert sd
            double deleteBeforeStdDevSum = 0;
            for (int b = 0; b < 5; ++b) {
                for (int q = 0; q < 63; ++q) {
                    double x = deleteBeforeMean - q;
                    deleteBeforeStdDevSum += result.qualityCountsBeforeDeletesPerBasePerPosition[b][q][position] * x * x;
                }
            }
            double deleteBeforeStdDev = std::sqrt(1.0 / deleteBeforeCount * deleteBeforeStdDevSum);
            if (deleteBeforeStdDevSum == 0)
                deleteBeforeStdDev = 0;
            printf("       %7.2f          %7.2f", deleteBeforeMean, deleteBeforeStdDev);
        }
        if (deleteAfterCount == 0) {
            printf("      %7s         %7s", "-", "-");
        } else {
            deleteAfterMean = deleteAfterSum / deleteAfterCount;
            // insert sd
            double deleteAfterStdDevSum = 0;
            for (int b = 0; b < 5; ++b) {
                for (int q = 0; q < 63; ++q) {
                    double x = deleteAfterMean - q;
                    deleteAfterStdDevSum += result.qualityCountsAfterDeletesPerBasePerPosition[b][q][position] * x * x;
                }
            }
            double deleteAfterStdDev = std::sqrt(1.0 / deleteAfterCount * deleteAfterStdDevSum);
            if (deleteAfterStdDevSum == 0)
                deleteAfterStdDev = 0;
            printf("      %7.2f         %7.2f", deleteAfterMean, deleteAfterStdDev);
        }
        // match/mismatch mean
        double mismatchSum = 0;
        double mismatchCount = 0;
        double mismatchMean = 0;
        double mismatchStdDev = 0;
        double matchSum = 0;
        double matchCount = 0;
        double matchMean = 0;
        double matchStdDev = 0;
        for (int b = 0; b < 4; ++b) {  // < 4 instead of < 5, ignoring matches/mismatches with N
            for (int a = 0; a < 4; ++a) {  // < 4 instead of < 5, ignoring matches/mismatches with N
                for (int q = 0; q < 63; ++q) {
                    if (a != b) {
                        mismatchSum += result.qualityCountsForMismatchPerMismatchPerPosition[a * 5 + b][q][position] * q;
                        mismatchCount += result.qualityCountsForMismatchPerMismatchPerPosition[a * 5 + b][q][position];
                    } else {
                        matchSum += result.qualityCountsForMismatchPerMismatchPerPosition[a * 5 + b][q][position] * q;
                        matchCount += result.qualityCountsForMismatchPerMismatchPerPosition[a * 5 + b][q][position];
                    }
                }
            }
        }
        if (mismatchCount > 0)
            mismatchMean = mismatchSum / mismatchCount;
        if (matchCount > 0)
            matchMean = matchSum / matchCount;
        if (mismatchCount > 0 || matchCount > 0) {
            // match/mismatch sd
            double mismatchStdDevSum = 0;
            double matchStdDevSum = 0;
            for (int b = 0; b < 4; ++b) {  // < 4 instead of < 5, ignoring matches/mismatches with N
                for (int a = 0; a < 4; ++a) {  // < 4 instead of < 5, ignoring matches/mismatches with N
                    for (int q = 0; q < 63; ++q) {
                        if (a != b) {
                            double x = mismatchMean - q;
                            mismatchStdDevSum += result.qualityCountsForMismatchPerMismatchPerPosition[a * 5 + b][q][position] * x * x;
                        } else {
                            double x = matchMean - q;
                            matchStdDevSum += result.qualityCountsForMismatchPerMismatchPerPosition[a * 5 + b][q][position] * x * x;
                        }
                    }
                }
            }
            if (mismatchCount > 0)
                mismatchStdDev = std::sqrt(1.0 / mismatchCount * mismatchStdDevSum);
            if (matchCount > 0)
                matchStdDev = std::sqrt(1.0 / matchCount * matchStdDevSum);
            if (mismatchCount == 0 && matchCount == 0)
                printf("  %7s       %7s %7s  %7s", "-", "-", "-", "-");
            else if (mismatchCount == 0 && matchCount > 0)
                printf("  %7s       %7s %7.2f  %7.2f", "-", "-", matchMean, matchStdDev);
            else if (mismatchCount > 0 && matchCount == 0)
                printf("  %7.2f       %7.2f %7s  %7s", mismatchMean, mismatchStdDev, "-", "-");
            else
                printf("  %7.2f       %7.2f %7.2f  %7.2f", mismatchMean, mismatchStdDev, matchMean, matchStdDev);
        }
        // total mean
        double totalSum = 0;
        double totalCount = 0;
        double totalMean = 0;
        for (int b = 0; b < 5; ++b) {
            for (int q = 0; q < 63; ++q) {
                totalSum += result.qualityCountsForInsertPerBasePerPosition[b][q][position] * q;
                totalCount += result.qualityCountsForInsertPerBasePerPosition[b][q][position];
                for (int a = 0; a < 5; ++a) {
                    totalSum += result.qualityCountsForMismatchPerMismatchPerPosition[a * 5 + b][q][position] * q;
                    totalCount += result.qualityCountsForMismatchPerMismatchPerPosition[a * 5 + b][q][position];
                }
            }
        }
        if (totalCount == 0) {
            printf("%7s   %7s", "-", "-");
        } else {
            totalMean = totalSum / totalCount;
            // total sd
            double totalStdDevSum = 0;
            for (int b = 0; b < 5; ++b) {
                for (int q = 0; q < 63; ++q) {
                    double x = totalMean - q;
                    totalStdDevSum += result.qualityCountsForInsertPerBasePerPosition[b][q][position] * x * x;
                    for (int a = 0; a < 5; ++a)
                        totalStdDevSum += result.qualityCountsForMismatchPerMismatchPerPosition[a * 5 + b][q][position] * x * x;
                }
            }
            double totalStdDev = std::sqrt(1 / totalCount * totalStdDevSum);
            printf("%7.2f   %7.2f", totalMean, totalStdDev);
        }
        std::cout << std::endl;

//         double total = length(fragmentStore.alignedReadStore);
//         printf("%9u   %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n", position, insertMean, insertStdDev, mismatchMean, mismatchSd, anyMean, anyStdDev, matchMean, matchStdDev, totalMean, totalStdDev);
    }
}

#endif  // #ifndef READ_ANALYZER_ANALYZE_ILLUMINA_H_
