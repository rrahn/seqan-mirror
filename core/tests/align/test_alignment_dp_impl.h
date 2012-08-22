// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_IMPL_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_IMPL_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/score.h>

#include <seqan/align/alignment_base.h>
#include <seqan/align/alignment_dp_band.h>
#include <seqan/align/alignment_dp_value.h>
#include <seqan/align/alignment_dp_formula.h>
#include <seqan/align/alignment_dp_matrix.h>
#include <seqan/align/alignment_dp_manager.h>
#include <seqan/align/alignment_dp_tracker.h>
#include <seqan/align/alignment_dp_impl.h>

#include <seqan/align/alignment_traceback.h>

using namespace seqan;


void testAlignmentDPImplBeginBandFirstPhase()
{
    DnaString seq = "ACGATCGACTACGTACGATCGACT"; //length = 24
    DnaString seq2 = "ACGTACGATCGACT"; //length = 14
    typedef Iterator<DnaString>::Type TIterator;

    DnaIterator itBegin = begin(seq);

    {  // Test for normal band.
        typedef Band<BandSwitchedOn<> > TBand;
        TIterator it = _beginBandFirstPhase(seq, seq2, TBand(-2,3));
        SEQAN_ASSERT(*it == 'A');
        SEQAN_ASSERT(it == itBegin);

        it = _beginBandFirstPhase(seq, seq2, TBand(0,3));
        SEQAN_ASSERT(*it == 'A');
        SEQAN_ASSERT(it == itBegin);

        it = _beginBandFirstPhase(seq, seq2, TBand(2,10));
        SEQAN_ASSERT(*it == 'G');
        SEQAN_ASSERT(it == itBegin + 2);
    }

    {  // Test for wide band.
        typedef Band<BandSwitchedOn<> > TBand;
        TIterator it = _beginBandFirstPhase(seq, seq2, TBand(-10,6));
        SEQAN_ASSERT(*it == 'A');
        SEQAN_ASSERT(it == itBegin);

        it = _beginBandFirstPhase(seq, seq2, TBand(0,18));
        SEQAN_ASSERT(*it == 'A');
        SEQAN_ASSERT(it == itBegin);

        it = _beginBandFirstPhase(seq, seq2, TBand(2,18));
        SEQAN_ASSERT(*it == 'G');
        SEQAN_ASSERT(it == itBegin + 2);
    }
}

void testAlignmentDPImplBeginBandMiddlePhase()
{

    DnaString seq = "ACGATCGACTACGTACGATCGACT";
    DnaString seq2 = "ACGTACGATCGACT";

    typedef Iterator<DnaString>::Type TIterator;

    {  // Test for normal band.
        typedef Band<BandSwitchedOn<> > TBand;
        TIterator it = _beginBandMiddlePhase(seq, seq2, TBand(-2,3));
        SEQAN_ASSERT(*it == 'A');
        SEQAN_ASSERT(it == begin(seq) + 3);

        it = _beginBandMiddlePhase(seq, seq2, TBand(0,10));
        SEQAN_ASSERT(*it == 'A');
        SEQAN_ASSERT(it == begin(seq) + 10);

        it = _beginBandMiddlePhase(seq, seq2, TBand(2,10));
        SEQAN_ASSERT(*it == 'A');
        SEQAN_ASSERT(it == begin(seq) + 10);
    }

    {  // Test for wide band.
        typedef Band<BandSwitchedOn<> > TBand;
        TIterator it = _beginWideBandMiddlePhase(seq, seq2, TBand(-10,6));
        SEQAN_ASSERT(*it == 'T');
        SEQAN_ASSERT(it == begin(seq) + 4);

        it = _beginWideBandMiddlePhase(seq, seq2, TBand(0,18));
        SEQAN_ASSERT(*it == 'A');
        SEQAN_ASSERT(it == begin(seq) + 14);

        it = _beginWideBandMiddlePhase(seq, seq2, TBand(2,18));
        SEQAN_ASSERT(*it == 'G');
        SEQAN_ASSERT(it == begin(seq) + 16);
    }
}

void testAlignmentDPImplBeginBandLastPhase()
{
    DnaString seq = "ACGATCGACTACGTACGATCGACT";
    DnaString seq2 = "ACGTACGATCGACT";

    typedef Iterator<DnaString>::Type TIterator;

    {  // Test for normal band.
        typedef Band<BandSwitchedOn<> > TBand;
        TIterator it = _beginBandEndPhase(seq, seq2, TBand(-2,3));
        SEQAN_ASSERT(*it == 'G');
        SEQAN_ASSERT(it == begin(seq) + 12);

        it = _beginBandEndPhase(seq, seq2, TBand(0,10));
        SEQAN_ASSERT(*it == 'A');
        SEQAN_ASSERT(it == begin(seq) + 14);

        it = _beginBandEndPhase(seq, seq2, TBand(2,10));
        SEQAN_ASSERT(*it == 'G');
        SEQAN_ASSERT(it == begin(seq) + 16);
    }

    {  // Test for wide band.
        typedef Band<BandSwitchedOn<> > TBand;
        TIterator it = _beginWideBandEndPhase(seq, seq2, TBand(-10,6));
        SEQAN_ASSERT(*it == 'G');
        SEQAN_ASSERT(it == begin(seq) + 6);

        it = _beginWideBandEndPhase(seq, seq2, TBand(0,18));
        SEQAN_ASSERT(*it == 'T');
        SEQAN_ASSERT(it == begin(seq) + 18);

        it = _beginWideBandEndPhase(seq, seq2, TBand(2,18));
        SEQAN_ASSERT(*it == 'T');
        SEQAN_ASSERT(it == begin(seq) + 18);
    }
}

void testAlignmentDPImplEndBand()
{
    DnaString seq = "ACGATCGACTACGTACGATCGACT";
    DnaString seq2 = "ACGTACGATCGACT";

    typedef Iterator<DnaString>::Type TIterator;

    {  // Test for normal band.
        typedef Band<BandSwitchedOn<> > TBand;
        TIterator it = _endBand(seq, seq2, TBand(-2,3));
        SEQAN_ASSERT(*it == 'A');
        SEQAN_ASSERT(it == begin(seq) + 17);

        it = _endBand(seq, seq2, TBand(0,10));
        SEQAN_ASSERT(it == end(seq));

        it = _endBand(seq, seq2, TBand(2,10));
        SEQAN_ASSERT(it == end(seq));
    }

    {  // Test for wide band.
        typedef Band<BandSwitchedOn<> > TBand;
        TIterator it = _endWideBand(seq, seq2, TBand(-10,6));
        SEQAN_ASSERT(*it == 'G');
        SEQAN_ASSERT(it == begin(seq) + 20);

        it = _endWideBand(seq, seq2, TBand(0,18));
        SEQAN_ASSERT(it == end(seq));

        it = _endWideBand(seq, seq2, TBand(2,18));
        SEQAN_ASSERT(it == end(seq));
    }
}



void testAlignmentDPImplIncrement()
{
    String<char> traceMatrix;
    resize(traceMatrix, 10, 1, Exact());
    typedef Iterator<String<char> >::Type TIterator;
    TIterator itBegin = begin(traceMatrix);
    TIterator it(itBegin);

    {
        increment(it, False());
        SEQAN_ASSERT(it == itBegin);
    }

    {
        increment(it, True());
        SEQAN_ASSERT(it == ++itBegin);
    }
}


void testAlignmentDPImplTrace()
{
    CharString traceString = "12341512513";
    CharIterator itBegin = begin(traceString);
    CharIterator it(itBegin);

    typedef DPValue<int, LinearGaps> TDPValue;
    TDPValue score = 10;
    {
        Tracker<int, CharIterator, ScoreOnly > traceTracker;
        track(traceTracker, score, it, True());

        SEQAN_ASSERT(traceTracker._maxScore == 10);
        SEQAN_ASSERT(length(traceTracker._tracePoints) == 0);

        score = 20;
        track(traceTracker, score, it, False());

        SEQAN_ASSERT(traceTracker._maxScore == 10);
        SEQAN_ASSERT(length(traceTracker._tracePoints) == 0);
    }
}

void testAlign2DPImplComputeCellLinear()
{
    typedef DPValue<int, LinearGaps> TDPValue;
    typedef Value<TDPValue>::Type TScoreValue;
    typedef Iterator<String<TraceBitMask::Type> >::Type TTraceIter;
    typedef Band<BandSwitchedOff> TBand;
    typedef Tracker<TScoreValue, TTraceIter, Default> TTracker;

    String<TraceBitMask::Type> traceMat;
    resize(traceMat, 1);
    traceMat[0] = +TraceBitMask::NONE;

    TTraceIter traceIter = begin(traceMat);
    DPFormula<Score<TScoreValue, Simple>, Global<> > recFormula(Score<TScoreValue, Simple>(2,-2,-4));

    TDPValue acticeCell;
    TDPValue prevD;
    TDPValue prevV;
    TDPValue prevH;

    setScore(prevD, 10);
    setScore(prevV, 12);
    setScore(prevH, 12);

    {
        TTracker traceTracker;
        typedef Cell<DPDirectionAll, True> TCell;
        computeCell(traceTracker, traceIter, acticeCell, prevD, prevV, prevH, 'A', 'C', recFormula, TCell(), True());
        SEQAN_ASSERT_EQ(traceIter, begin(traceMat));
        SEQAN_ASSERT_EQ(*traceIter, (+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL | +TraceBitMask::VERTICAL));
        SEQAN_ASSERT_EQ(traceTracker._maxScore, 8);
        SEQAN_ASSERT_EQ(value(traceTracker._tracePoints[0]), (+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL | +TraceBitMask::VERTICAL));
    }


    {
        *traceIter = +TraceBitMask::NONE;
        TTracker traceTracker;
        typedef Cell<DPDirectionAll, True> TCell;
        computeCell(traceTracker, traceIter, acticeCell, prevD, prevV, prevH, 'A', 'A', recFormula, TCell(), False());
        SEQAN_ASSERT_EQ(traceIter, begin(traceMat));
        SEQAN_ASSERT_EQ(*traceIter, (+TraceBitMask::NONE));
        SEQAN_ASSERT_EQ(traceTracker._maxScore, 12);
        SEQAN_ASSERT_EQ(value(traceTracker._tracePoints[0]), (+TraceBitMask::NONE));
    }
}

void testAlign2DPImplComputeCellAffine()
{
    typedef DPValue<int, AffineGaps> TDPValue;
    typedef Value<TDPValue>::Type TScoreValue;
    typedef Iterator<String<TraceBitMask::Type> >::Type TTraceIter;
    typedef Band<BandSwitchedOff> TBand;
    typedef Tracker<TScoreValue, TTraceIter, Default> TTracker;

    String<TraceBitMask::Type> traceMat;
    resize(traceMat, 1);
    traceMat[0] = +TraceBitMask::NONE;

    TTraceIter traceIter = begin(traceMat);
    DPFormula<Score<TScoreValue, Simple>, Global<> > recFormula(Score<TScoreValue, Simple>(2,-2,-4));

    TDPValue acticeCell;
    TDPValue prevD;
    TDPValue prevV;
    TDPValue prevH;

    setScore(prevD, 10);
    setScore(prevV, 12);
    setScore(prevH, 12);

    {
        TTracker traceTracker;
        typedef Cell<DPDirectionAll, True> TCell;
        computeCell(traceTracker, traceIter, acticeCell, prevD, prevV, prevH, 'A', 'C', recFormula, TCell(), True());
        SEQAN_ASSERT_EQ(traceIter, begin(traceMat));
        SEQAN_ASSERT_EQ(*traceIter, (+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL_OPEN | +TraceBitMask::VERTICAL_OPEN));
        SEQAN_ASSERT_EQ(traceTracker._maxScore, 8);
        SEQAN_ASSERT_EQ(value(traceTracker._tracePoints[0]), (+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL_OPEN | +TraceBitMask::VERTICAL_OPEN));
    }


    {
        *traceIter = +TraceBitMask::NONE;
        TTracker traceTracker;
        typedef Cell<DPDirectionAll, True> TCell;
        computeCell(traceTracker, traceIter, acticeCell, prevD, prevV, prevH, 'A', 'A', recFormula, TCell(), False());
        SEQAN_ASSERT_EQ(traceIter, begin(traceMat));
        SEQAN_ASSERT_EQ(*traceIter, (+TraceBitMask::NONE));
        SEQAN_ASSERT_EQ(traceTracker._maxScore, 12);
        SEQAN_ASSERT_EQ(value(traceTracker._tracePoints[0]), (+TraceBitMask::NONE));
    }
}

template <typename TDPValue>
void testAlign2DPImplFillColumnUnbanded()
{
    typedef typename Value<TDPValue>::Type TScoreValue;
    typedef typename Iterator<String<typename TraceBitMask::Type> >::Type TTraceIter;
    typedef Band<BandSwitchedOff> TBand;
    typedef typename Iterator<DPMatrix<TDPValue> >::Type TDPIter;
    typedef typename Iterator<CharString>::Type TCharIterator;
    typedef AlignmentProfile<Global<>, typename Spec<TDPValue>::Type, Traceback<Default> > TAlignmentProfile;
    typedef DPManager<TAlignmentProfile, TBand> TDPManager;
    typedef typename IsTracebackOn<TAlignmentProfile>::Type TIsTraceback;

    typedef Tracker<TScoreValue, TTraceIter, typename GetTrackerSpec<TAlignmentProfile>::Type > TTracker;

    CharString seq0 = "Hello Test";
    CharString seq1 = "Hello Seqan";

    TBand band;
    DPMatrix<TDPValue> dpMatrix(seq0, seq1, band);

    String<typename TraceBitMask::Type> traceMat;
    resize(traceMat, (length(seq0)+ 1)*(length(seq1)+1));


    TDPManager colManager(getColumnSize(seq0, seq1, band), getInitialColumnSize(seq0, seq1, band));

    TDPIter activeIt = begin(dpMatrix, colManager);
    TDPIter endIt = end(dpMatrix) - 1;
    TTraceIter traceIter = begin(traceMat, colManager);
    DPFormula<Score<TScoreValue, Simple>, Global<> > recFormula(Score<TScoreValue, Simple>(2,-2,-2));
    TCharIterator seqVIter = begin(seq1) + colManager._spanSeqVBegin;
    TCharIterator seqVStop = begin(seq1) + colManager._spanSeqVEnd;
    TTracker traceTracker;

    { // InitialColumn
        update(colManager, DPNoBandInitPhase());
        fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0 [0], recFormula,
                   DPNoBandInitPhase(), colManager);
        SEQAN_ASSERT_EQ(getScore(*activeIt), -22);
        SEQAN_ASSERT(activeIt == endIt);
    }

    seqVIter = begin(seq1) + colManager._spanSeqVBegin;
    seqVStop = begin(seq1) + colManager._spanSeqVEnd;
    activeIt += colManager._spanDp;
    traceIter += colManager._spanTrace;
    { // SeconDPPhase
        update(colManager, DPNoBandPhase());
        fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[0], recFormula, DPNoBandPhase(), TDPManager());
        SEQAN_ASSERT_EQ(getScore(*activeIt), -18);
        SEQAN_ASSERT(activeIt == endIt);

    }

    seqVIter = begin(seq1) + colManager._spanSeqVBegin;
    seqVStop = begin(seq1) + colManager._spanSeqVEnd;
    activeIt += colManager._spanDp;
    traceIter += colManager._spanTrace;
    { // SeconDPPhase
        update(colManager, DPNoBandPhase());
        fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[1], recFormula, DPNoBandPhase(), TDPManager());
        SEQAN_ASSERT_EQ(getScore(*activeIt), -14);
        SEQAN_ASSERT(activeIt == endIt);
    }
}

template <typename TColIt>
void printColumn(TColIt const & itEnd, TColIt const & begin)
{
    TColIt it = begin;
    for (; it != itEnd; ++it)
    {
//        std::cout << getScore(*it) << "\t";
    }
//    std::cout << getScore(*it) << "\n";
}

template <typename TTraceValue>
void printTraceMatrix(String<TTraceValue> const & traceMat)
{
    typedef typename Iterator<String<TTraceValue> >::Type TIter;

    int colSize = 7;
    int rowSize = 11;


    for (int row = 0; row < colSize; ++row)
    {
        for (int col = 0; col < rowSize; ++col)
        {
            std::cout << (int) traceMat[col*colSize + row] << "\t";
        }
        std::cout << "\n";
    }

}


template <typename TDPValue>
void testAlign2DPImplFillColumnBanded()
{
    typedef typename Value<TDPValue>::Type TScoreValue;
    typedef typename Iterator<String<typename TraceBitMask::Type> >::Type TTraceIter;
    typedef Band<BandSwitchedOn<> > TBand;
    typedef typename Iterator<DPMatrix<TDPValue> >::Type TDPIter;
    typedef Tracker<TScoreValue, TTraceIter, Default > TTracker;
    typedef typename Iterator<CharString>::Type TCharIterator;
    typedef DPManager<AlignmentProfile<Global<>, typename Spec<TDPValue>::Type, Traceback<Default> >,  TBand> TDPManager;


    CharString seq0 = "HELLO TEST";
    CharString seq1 = "HELLO SEQAN";

    TBand band(-3,3);
    DPMatrix<TDPValue> dpMatrix(seq0, seq1, band);

    String<typename TraceBitMask::Type> traceMat;
    resize(traceMat, (length(seq0)+1) * 7);

    DPFormula<Score<TScoreValue, Simple>, Global<> > recFormula(Score<TScoreValue, Simple>(2,-2,-2));
    size_t colSize = getColumnSize(seq0, seq1, band);
    size_t initialColSize = getInitialColumnSize(seq0, seq1, band);
    TDPManager columnManager(colSize, initialColSize);

    TCharIterator seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    TCharIterator seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    TDPIter activeIt = begin(dpMatrix, columnManager);
    TDPIter endIt = end(dpMatrix) - 1;
    TTraceIter traceIter = begin(traceMat, columnManager);
    TTracker traceTracker;
    TDPIter beginMatrix = begin(dpMatrix);
    // InitialColumn
    update(columnManager, DPBandInitPhase());
    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[0], recFormula, DPBandInitPhase(), columnManager);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), -6);
    SEQAN_ASSERT(activeIt == endIt);
    printColumn(activeIt, beginMatrix);

    // BandOpenColumn
    update(columnManager, DPBandFirstPhase());
    seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    activeIt += columnManager._spanDp;
    traceIter += columnManager._spanTrace;

    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[0], recFormula, DPBandFirstPhase(), columnManager);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), -4);
    SEQAN_ASSERT(activeIt == endIt);

    printColumn(activeIt, beginMatrix);

    update(columnManager, DPBandFirstPhase());
    seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    activeIt += columnManager._spanDp;
    traceIter += columnManager._spanTrace;

    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[1], recFormula, DPBandFirstPhase(), columnManager);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), -2);
    SEQAN_ASSERT(activeIt == endIt);

    printColumn(activeIt, beginMatrix);

    update(columnManager, DPBandFirstPhase());
    seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    activeIt += columnManager._spanDp;
    traceIter += columnManager._spanTrace;

    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[2], recFormula, DPBandFirstPhase(), columnManager);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), 0);
    SEQAN_ASSERT(activeIt == endIt);

    printColumn(activeIt, beginMatrix);

    // Full Column
    update(columnManager, DPBandMiddlePhase());
    seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    activeIt += columnManager._spanDp;
    traceIter += columnManager._spanTrace;

    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[3], recFormula, DPBandMiddlePhase(), columnManager);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), 2);
    SEQAN_ASSERT(activeIt == endIt);

    printColumn(activeIt, beginMatrix);

    update(columnManager, DPBandMiddlePhase());
    seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    activeIt += columnManager._spanDp;
    traceIter += columnManager._spanTrace;

    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[4], recFormula, DPBandMiddlePhase(), columnManager);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), 4);
    SEQAN_ASSERT(activeIt == endIt);

    printColumn(activeIt, beginMatrix);

    update(columnManager, DPBandMiddlePhase());
    seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    activeIt += columnManager._spanDp;
    traceIter += columnManager._spanTrace;

    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[5], recFormula, DPBandMiddlePhase(), columnManager);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), 6);
    SEQAN_ASSERT(activeIt == endIt);

    printColumn(activeIt, beginMatrix);

    update(columnManager, DPBandMiddlePhase());
    seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    activeIt += columnManager._spanDp;
    traceIter += columnManager._spanTrace;

    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[6], recFormula, DPBandMiddlePhase(), columnManager);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), 4);
    SEQAN_ASSERT(activeIt == endIt);

    printColumn(activeIt, beginMatrix);

    update(columnManager, DPBandMiddlePhase());
    seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    activeIt += columnManager._spanDp;
    traceIter += columnManager._spanTrace;

    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[7], recFormula, DPBandMiddlePhase(), columnManager);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), 6);
    SEQAN_ASSERT(activeIt == endIt);

    printColumn(activeIt, beginMatrix);

    // BandCloseColumn
    _correctSpan(columnManager);
    update(columnManager, DPBandLastPhase());
    seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    activeIt += columnManager._spanDp;
    traceIter += columnManager._spanTrace;

    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[8], recFormula, DPBandLastPhase(), columnManager);
    printColumn(activeIt, beginMatrix);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), 6);
    SEQAN_ASSERT(activeIt == --endIt);


    update(columnManager, DPBandLastPhase());
    seqVIter = begin(seq1) + columnManager._spanSeqVBegin;
    seqVStop = begin(seq1) + columnManager._spanSeqVEnd;
    activeIt += columnManager._spanDp;
    traceIter += columnManager._spanTrace;

    fillColumn(traceTracker, traceIter, activeIt, seqVIter, seqVStop, seq0[9], recFormula, DPBandLastPhase(), columnManager);
    SEQAN_ASSERT_EQ(getScore(*(activeIt)), 6);
    SEQAN_ASSERT(activeIt == --endIt);

    printColumn(activeIt, beginMatrix);

}

void testAlign2DPImplIsValidDPSettings()
{
    typedef AlignmentProfile<Global<>, LinearGaps, TracebackSwitchedOff > TAlignProfile1;
    typedef AlignmentProfile<Global<FreeEndGaps<True, True, True, True> >, AffineGaps, TracebackSwitchedOn> TAlignProfile2;
    typedef AlignmentProfile<Local<>, AffineGaps, TracebackSwitchedOn> TAlignProfile3;

    {  // Check result when using empty sequences.
        DnaString seqH = "";
        Dna5String seqV = "ACGTGACT";

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile2()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile3()), false);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-3,3), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-3,3), TAlignProfile2()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-3,3), TAlignProfile3()), false);

        seqH = "ACGTACGTAC";
        resize(seqV, 0);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile2()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile3()), false);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-3,3), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-3,3), TAlignProfile2()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-3,3), TAlignProfile3()), false);

        resize(seqH, 0);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile2()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile3()), false);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-3,3), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-3,3), TAlignProfile2()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-3,3), TAlignProfile3()), false);
    }

    {  // Check empty intersection between band and DP matrix.
        DnaString seqH =  "ACGTGACT";
        Dna5String seqV = "ACGTGACT";

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile1()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile3()), true);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-10,-9), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-10,-9), TAlignProfile2()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-10,-9), TAlignProfile3()), false);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(9,10), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(9,10), TAlignProfile2()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(9,10), TAlignProfile3()), false);
    }

    {  // Check band begins before beginning of horizontal sequence.
        DnaString seqH =  "ACGTGACT";
        Dna5String seqV = "ACGTGACT";

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile1()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile3()), true);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-4, 0), TAlignProfile1()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-4, 0), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-4, 0), TAlignProfile3()), true);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-10,-1), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-10,-1), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-10,-1), TAlignProfile3()), true);
    }

    {  // Check band begins before beginning of vertical sequence.
        DnaString seqH =  "ACGTGACT";
        Dna5String seqV = "ACGTGACT";

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile1()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile3()), true);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(0,4), TAlignProfile1()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(0,4), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(0,4), TAlignProfile3()), true);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(1,10), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(1,10), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(1,10), TAlignProfile3()), true);
    }

    {  // Check band ends behind end of vertical sequence.
        DnaString seqH =  "ACGTGACTACGTGACT";
        Dna5String seqV = "ACGTGACT";

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile1()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile3()), true);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-4, 8), TAlignProfile1()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-4, 8), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-4, 8), TAlignProfile3()), true);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-10, 7), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-10, 7), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-10, 7), TAlignProfile3()), true);
    }

    {  // Check band ends behind end of horizontal sequence.
        DnaString seqH =  "ACGTGACT";
        Dna5String seqV = "ACGTGACTACGTGACT";

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile1()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOff>(), TAlignProfile3()), true);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-8,4), TAlignProfile1()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-8,4), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-8,4), TAlignProfile3()), true);

        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-7,10), TAlignProfile1()), false);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-7,10), TAlignProfile2()), true);
        SEQAN_ASSERT_EQ( _isValidDPSettings(seqH, seqV, Band<BandSwitchedOn<> >(-7,10), TAlignProfile3()), true);
    }


}

template <typename TAlgoTag, typename TSeqH, typename TSeqV>
int testAlignmentDPImplBand(int lowerDiagonal, int upperDiagonal, TSeqH const & seqH, TSeqV const & seqV)
{
    typedef TraceSegment<size_t, size_t> TTraceSegment;
    typedef Band<BandSwitchedOn<> > TBand;
    typedef Iterator<String<char> >::Type TIter;
    typedef AlignmentProfile<TAlgoTag, LinearGaps, TracebackSwitchedOff> TAlignmentProfile;

    String<TTraceSegment> trace;
    Score<int, Simple> score(2, -2, -2);
    Tracker<int, TIter, Default> tracker;
    TBand band(lowerDiagonal, upperDiagonal);

    CharString traceMatrix;
    resize(traceMatrix, getTraceMatrixSize(seqH, seqV, band, False()));

//    std::cout << seqH << std::endl;

    computeDPMatrix(tracker, traceMatrix, seqH, seqV, score, band, TAlignmentProfile());
//    std::cout << "\n" << std::endl;
//    std::cout << tracker._maxScore << std::endl;
    return tracker._maxScore;
}

// Testing bands where upper diagonal crosses the horizontal sequence and the lower diagonal
// crosses the vertical sequence.
SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_noright_end_nottop_noleft_bottom_right)
{
    // small band and wide band
    //x x x x
    //x       x
    //x         x
    //  x         x
    //    x         x
    //      x       x
    //        x x x x

    Dna5String seqH = "AAAACCCCGGGGTTTT";
    DnaString seqV  = "AAAAGGGGCCCCTTTT";

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,1, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,2, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-2,1, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,3, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-3,1, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-4,1, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,4, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-4,4, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-6,6, seqH, seqV), 8);


    // wide band
    //x x x x x
    //x         x
    //x           x
    //x             x
    //  x             x
    //    x           x
    //      x x x x x x

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-7,7, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-8,3, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-3,8, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-11,11, seqH, seqV), 8);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_noright_end_nottop_noleft_bottom_noright)
{
    // small band and wide band
    //x x x x
    //x       x
    //x         x
    //  x         x
    //    x         x
    //      x         x
    //        x x x x x x

    Dna5String seqH = "AAAAGGGGCCCCTTTT";
    DnaString seqV  = "AAAACCCC";
    typedef FreeEndGaps<False, False, True, False> TFreeEndGaps;
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,1, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-1,1, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-2,1, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-5,2, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-6,3, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-7,1, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-1,8, seqH, seqV), 8);

    // wide band
    //x x x x x
    //x         x
    //x           x
    //x             x
    //  x             x
    //    x             x
    //      x x x x x x x x

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-7,6, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-7,2, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-6,3, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-5,5, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-7,8, seqH, seqV), 8);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_noright_end_nottop_noleft_nobottom_right)
{
    // small band and wide band
    //x x x x
    //x       x
    //x         x
    //  x       x
    //    x     x
    //      x   x
    //        x x
    //          x
    DnaString seqH  = "AAAACCCC";
    Dna5String seqV = "AAAAGGGGCCCCTTTT";
    typedef FreeEndGaps<False, False, False, True> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,1, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-1,2, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-2,1, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-3,2, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-4,3, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-1,7, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-2,7, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-3,7, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-4,7, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-8,7, seqH, seqV), 8);
}

// Testing bands where upper diagonal crosses the horizontal sequence and the lower diagonal
// crosses the horizontal sequence.
SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_noleft_nobottom_noright_end_nottop_noleft_bottom_right)
{
    // small band and wide band
    //x x x x x x
    //  x         x
    //    x         x
    //      x         x
    //        x         x
    //          x       x
    //            x x x x

    DnaString seqH  = "AAAACCCCGGGGTTTT";
    Dna5String seqV = "AAAAGGGG";

    typedef FreeEndGaps<True, False, True, True> TFreeEndGaps;
    typedef FreeEndGaps<True, False, False, False> TFreeEndGaps2;
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps2> >(1,9, seqH, seqV), -2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(1,9, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(1,9, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(7,9, seqH, seqV), -12);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps2> >(5,10, seqH, seqV), -10);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(5,10, seqH, seqV), -4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(7,15, seqH, seqV), -2);

    // wide band
    //x x x x x x x x
    //  x             x
    //    x             x
    //      x             x
    //        x             x
    //          x           x
    //            x x x x x x

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,9, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(5,15, seqH, seqV), -2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,15, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(1,15, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(2,15, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(3,15, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(4,15, seqH, seqV), 0);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_noleft_nobottom_noright_end_nottop_noleft_bottom_noright)
{
    // small band and wide band
    //x x x x x x
    //  x         x
    //    x         x
    //      x         x
    //        x         x
    //          x         x
    //            x x x x x x

    DnaString seqH  = "AAAACCCCGGGGTTTTGGGG";
    Dna5String seqV = "AAAAGGGG";
    typedef FreeEndGaps<True, False, True, False> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,1, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,1, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,2, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,3, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,4, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,8, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(8,12, seqH, seqV), 0);

    // wide band
    //x x x x x x x x
    //  x             x
    //    x             x
    //      x             x
    //        x             x
    //          x             x
    //            x x x x x x x x

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(0,12, seqH, seqV), -8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(1,12, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,9, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(3,12, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,12, seqH, seqV), 8);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_noleft_nobottom_noright_end_nottop_noleft_nobottom_right)
{
    // small band and wide band
    //x x x x x x
    //  x         x
    //    x         x
    //      x       x
    //        x     x
    //          x   x
    //            x x
    //              x

    DnaString seqH  = "AAAACCCCGGGG";
    Dna5String seqV = "AAAAGGGG";

    typedef FreeEndGaps<True, False, False, True> TFreeEndGaps;


    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(5,10, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(4,5, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(4,11, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(5,11, seqH, seqV), -2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(6,11, seqH, seqV), -2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(10,11, seqH, seqV), -2);
}

// Testing bands where upper diagonal crosses the horizontal sequence somewhere behind the
// end of the horizontal sequence and the lower diagonal crosses the horizontal sequence.
SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_noleft_nobottom_right_end_top_noleft_bottom_right)
{
    // only wide band
    //x x x x x x x x x x
    //  x               x
    //    x             x
    //      x           x
    //        x         x
    //          x       x
    //            x x x x

    DnaString seqH  = "AAAACCCCGGGGTTTT";
    Dna5String seqV = "AAAAGGGG";

    typedef FreeEndGaps<True, False, True, True> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(1,16, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(0,16, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,16, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(1,16, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(2,16, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(3,16, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(4,16, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,17, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(1,17, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(2,17, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(3,17, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(4,17, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,1000, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(1,1000, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(2,1000, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(3,1000, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(4,1000, seqH, seqV), 0);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_noleft_nobottom_right_end_top_noleft_nobottom_right)
{
    // small band and wide band (only FirstRow needed)
    //x x x x x x
    //  x       x
    //    x     x
    //      x   x
    //        x x
    //          x

    DnaString seqH  = "AAAACCCCGGGGTTTT";
    Dna5String seqV = "AAAATTTT";

    typedef FreeEndGaps<True, False, False, True> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(8,16, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(8,16, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(9,16, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(15,16, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(8,17, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(9,17, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(8,1000, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(9,1000, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(16,17, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(16,1000, seqH, seqV), 0);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_right_end_top_noleft_bottom_right)
{
    // only wide band
    // x x x x x x x x
    // x             x
    // x             x
    //   x           x
    //     x         x
    //       x       x
    //         x x x x

    Dna5String seqH = "AAAACCCCGGGGTTTT";
    DnaString seqV  = "AAAAGGGGCCCCTTTT";

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,16, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,17, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,1000, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-2,16, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-2,17, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-2,1000, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-15,16, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-15,17, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-15,1000, seqH, seqV), 8);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_right_end_top_noleft_nobottom_right)
{
    // small band and wide band (only FirstRow needed)
    //  x x x x
    //  x     x
    //  x     x
    //    x   x
    //      x x
    //        x

    Dna5String seqH = "AAAATTTT";
    DnaString seqV  = "AAAAGGGGCCCCTTTT";

    typedef FreeEndGaps<False, False, False, True> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-7,8, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-8,8, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-7,9, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-7,1000, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-7,8, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-7,9, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-7,1000, seqH, seqV), 0);

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-1,8, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-1,9, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-1,1000, seqH, seqV), 0);
}

// Testing bands where upper diagonal crosses the vertical sequence and the lower diagonal
// crosses the vertical sequence.
SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_notop_left_nobottom_norigth_end_notop_noleft_bottom_right)
{
    // only small band
    // x
    // x x
    // x   x
    // x     x
    //   x     x
    //     x     x
    //       x x x

    Dna5String seqH = "CCCCTTTT";
    DnaString seqV  = "AAAACCCCGGGGTTTT";

    typedef FreeEndGaps<False, True, False, False> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-9,-1, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-15,-1, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-9,-1, seqH, seqV), 8);  // BUG
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-15,-1, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-9,-7, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-9,-6, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-9,-5, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-9,-4, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-15,-7, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-15,-6, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-15,-5, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-15,-4, seqH, seqV), 8);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_notop_left_nobottom_noright_end_notop_noleft_bottom_noright)
{
    // x
    // x x
    // x   x
    // x     x
    //   x     x
    //     x     x
    //       x x x x

    Dna5String seqH = "AAAACCCCGGGGTTTT";
    DnaString seqV  = "AAAAGGGGCCCCTTTT";

    typedef FreeEndGaps<False, True, True, False> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,0, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-4, -2, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-1,0, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-2,0, seqH, seqV), 4);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-3,0, seqH, seqV), 6);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-4,0, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-15, 0, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-15, -14, seqH, seqV), -2);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_notop_left_nobottom_noright_end_notop_noleft_nobottom_right)
{
    // x
    // x x
    // x   x
    // x     x
    //   x   x
    //     x x
    //       x

    Dna5String seqH = "AAAACCCC";
    DnaString seqV  = "AAAAGGGGCCCCTTTT";

    typedef FreeEndGaps<False, True, False, True> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-8,0, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-8, -1, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-8,0, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-8, -7, seqH, seqV), -12);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-1, 0, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-6, -2, seqH, seqV), 4);
}

// horizontal begin pos on left, vertical begin pos behind seqV
SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_notop_left_nobottom_noright_end_notop_left_bottom_noright)
{
    // x
    // x x
    // x   x
    // x     x
    // x       x
    // x         x
    // x x x x x x x

    Dna5String seqH = "AAAACCCCGGGGTTTT";
    DnaString seqV  = "AAAAGGGGCCCCTTTT";

    typedef FreeEndGaps<False, True, True, False> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-16,0, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-20000,0, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-16, -2, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-20000, -2, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-16,0, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000,0, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-16, -4, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000, -4, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-16, -11, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000, -11, seqH, seqV), 0);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_notop_left_nobottom_noright_end_notop_left_bottom_right)
{
    // x
    // x x
    // x   x
    // x     x
    // x     x
    // x     x
    // x x x x

    Dna5String seqH = "AAAACCCC";
    DnaString seqV  = "AAAAGGGGCCCCTTTT";

    typedef FreeEndGaps<False, True, False, False> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-16,0, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-20000,0, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-16, -2, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-20000, -2, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-16,0, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000,0, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-16, -1, seqH, seqV), -2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000, -1, seqH, seqV), -2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-16, -7, seqH, seqV), -14);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000, -7, seqH, seqV), -14);
}

// horizontal begin pos on top, vertical begin pos behind seqV
SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_noright_end_notop_left_bottom_noright)
{
    // only wide band
    // x x x
    // x     x
    // x       x
    // x         x
    // x           x
    // x             x
    // x x x x x x x x x

    Dna5String seqH = "AAAACCCCGGGGTTTT";
    DnaString seqV  = "AAAAGGGG";

    typedef FreeEndGaps<False, False, True, False> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-8,8, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-20000,8, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-8, 6, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-20000, 6, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-8,1, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000,1, seqH, seqV), 2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-8, 8, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000, 8, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-8, 4, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000, 4, seqH, seqV), 8);
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_noright_end_notop_left_bottom_right)
{
    // only wide band
    // x x x
    // x     x
    // x       x
    // x         x
    // x           x
    // x           x
    // x x x x x x x

    Dna5String seqH = "AAAACCCCGGGGTTTT";
    DnaString seqV  = "AAAAGGGG";

    typedef FreeEndGaps<False, False, True, False> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-8,9, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-20000,9, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-8,9, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000,9, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-8, 15, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000, 15, seqH, seqV), 8);
}

// horizontal begin pos behind seqH, vertical begin pos behind seqV
SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_top_left_bottom_right_end_top_left_bottom_right)
{
    // only wide band
    // x x x x x x x
    // x           x
    // x           x
    // x           x
    // x           x
    // x           x
    // x x x x x x x

    Dna5String seqH = "AAAACCCCGGGGTTTT";
    DnaString seqV  = "AAAAGGGG";

    typedef FreeEndGaps<False, False, True, False> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-8,16, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-20000,20000, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-8,16, seqH, seqV), 8);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000,20000, seqH, seqV), 8);
}

// horizontal begin pos behind seqH, vertical begin pos behind seqV
SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_begin_notop_noleft_nobottom_noright_end_notop_noleft_nobottom_noright)
{
    // no band ...
    Dna5String seqH = "AAAACCCCGGGGTTTT";
    DnaString seqV  = "AAAAGGGG";

    typedef FreeEndGaps<True, True, True, True> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-9,-9, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-20000,-9, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(17,17, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(17,2000, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-9,-9, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-20000,-9, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(17,17, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(17,2000, seqH, seqV), MinValue<int>::VALUE);

}

// horizontal begin pos = vertical begin pos
SEQAN_DEFINE_TEST(test_alignment_dp_impl_band_width_one_diagonal)
{
    // special case.
    // x
    //   x
    //     x
    //       x

    Dna5String seqH = "AAAACCCCGGGGTTTT";
    DnaString seqV  = "AAAAGGGG";

    typedef FreeEndGaps<True, True, True, True> TFreeEndGaps;

    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-8,-8, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(-1,-1, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(0,0, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(1,1, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(8,8, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<> >(16,16, seqH, seqV), MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-8,-8, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(-1,-1, seqH, seqV), -2);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(0,0, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(1,1, seqH, seqV), 0);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(8,8, seqH, seqV), -16);
    SEQAN_ASSERT_EQ(testAlignmentDPImplBand<Global<TFreeEndGaps> >(16,16, seqH, seqV), 0);
}


SEQAN_DEFINE_TEST(test_alignment_dp_impl_begin_first_phase)
{

    testAlignmentDPImplBeginBandFirstPhase();
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_begin_middle_phase)
{
    testAlignmentDPImplBeginBandMiddlePhase();
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_begin_last_phase)
{
    testAlignmentDPImplBeginBandLastPhase();
}


SEQAN_DEFINE_TEST(test_alignment_dp_impl_end_band)
{
    testAlignmentDPImplEndBand();
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_increment)
{
    testAlignmentDPImplIncrement();
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_trace)
{
    testAlignmentDPImplTrace();
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_compute_cell)
{
    testAlign2DPImplComputeCellLinear();
    testAlign2DPImplComputeCellAffine();
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_fill_column)
{
    testAlign2DPImplFillColumnUnbanded<DPValue<int, LinearGaps> >();
    testAlign2DPImplFillColumnUnbanded<DPValue<int, AffineGaps> >();
    testAlign2DPImplFillColumnBanded<DPValue<int, LinearGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_impl_is_valid_dp_settings)
{
    testAlign2DPImplIsValidDPSettings();
}

#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_IMPL_H_
