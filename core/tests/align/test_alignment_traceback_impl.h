// ==========================================================================
//                      test_alignment_traceback_impl.h
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

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_TRACEBACK_IMPL_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_TRACEBACK_IMPL_H_

#include <seqan/basic.h>
#include <seqan/score.h>

#include <seqan/align/alignment_dp_value.h>
#include <seqan/align/alignment_dp_band.h>
#include <seqan/align/alignment_dp_matrix.h>
#include <seqan/align/alignment_dp_formula.h>
#include <seqan/align/alignment_dp_tracker.h>
#include <seqan/align/alignment_dp_manager.h>
#include <seqan/align/alignment_dp_impl.h>
#include <seqan/align/alignment_dp.h>

#include <seqan/align/alignment_traceback_tracesegment.h>
#include <seqan/align/alignment_traceback_impl.h>
#include <seqan/align/alignment_traceback_adaptor.h>
#include <seqan/align/alignment_traceback.h>


using namespace seqan;

//template <typename TAlgorithm, typename TValue1>
//void testAlignmentTracebackImplTracebackComputeTraceback(TValue1 const & nominalVal1)
//{
//    //the trace iter is the iterator to the cell with the maximum value within the tracematrix.. (independent of local and global)
//    //the trace can be computed a following TraceOne, TraceAll - two specifier to record one or all traces.
//    //record the traceback ...
//        //file, Align Object, stream, AlingmentGraph ... segment wise (longest common subsequence) |
////    computeTraceback(traceIter, seq0, seq1, TAlgorithm());
//}

void testAlignmentTracebackImplTraceCoreGetTraceMatrixSize()
{

    DnaString str0 = "ACGT";
    DnaString str1 = "ACGTTGCA";

    SEQAN_ASSERT_EQ(getTraceMatrixSize(str0, str1, Band<BandSwitchedOff>(), True()), 45u);
    SEQAN_ASSERT_EQ(getTraceMatrixSize(str0, str1, Band<BandSwitchedOff>(), False()), 1u);
    SEQAN_ASSERT_EQ(getTraceMatrixSize(str0, str1, Band<BandSwitchedOn<> >(-2,3), False()), 1u);
    SEQAN_ASSERT_EQ(getTraceMatrixSize(str0, str1, Band<BandSwitchedOn<> >(-2,3), True()), 30u);
}


// SEQAN_DEFINE_TEST(test_align2_traceback_compute_alignments)
// {
//
//     typedef AlignmentProfile<Global<LinearGaps, FullMatrix, Traceback<SingleTrace> >, FreeEndGaps<> > TAlignmentProfile;
//
//     DnaString seq0 = "ACGGCGGCGCGATGACGACT";
//     DnaString seq1 = "ACGACGCACGTACGACGTACTCG";
//
//     Score<int, Simple> scoreScheme(2,-1,-4);
//
//     Align<DnaString> alignObj;
//     resize(rows(alignObj), 2);
//     assignSource(row(alignObj, 0), seq0);
//     assignSource(row(alignObj, 1), seq1);
//     StringSet<DnaString> str;
//     appendValue(str, seq0);
//     appendValue(str, seq1);
// //    StringSet<TraceSegment<size_t, size_t> > traceSegments;
// //    doAlign(traceSegments, seq0, seq1, scoreScheme, Band<Off>(), TAlignmentProfile());
// //
// //    Nothing target;
// //    adaptTraceSegments(target, seq0, seq1, traceSegments);
//
//
// //    _globalAlignment(TAlign& align,
// //                     TStringSet const& str,
// //                     TScore const& sc,
// //                     TAlignConfig const,
// //                     NeedlemanWunsch)
//
//     typedef typename Size<DnaString>::Type TSize;
//     AlignTraceback<TSize> trace;
//     std::cout << "Score: " << globalAlignment(trace, str, scoreScheme, AlignConfig<>(), NeedlemanWunsch()) << std::endl;
//     _pumpTraceToAlign(alignObj, trace);
//     std::cout << alignObj << std::endl;
//
//     clearGaps(alignObj);
//     clearClipping(alignObj);
//     AlignTraceback<TSize> trace2;
//
//     std::cout << "Score: " << globalAlignment(trace2, str, scoreScheme, -7, 7, BandedNeedlemanWunsch()) << std::endl;
//     _pumpTraceToAlign(alignObj, trace2);
//     std::cout << alignObj << std::endl;
//
//     clearGaps(alignObj);
//     clearClipping(alignObj);
//     AlignTraceback<TSize> trace3;
//     std::cout << "Score: " << globalAlignment(trace3, str, Score<int, Simple>(2, -4, -1, -2), AlignConfig<>(), Gotoh()) << std::endl;
//     _pumpTraceToAlign(alignObj, trace3);
//     std::cout << alignObj << std::endl;
//
//     clearGaps(alignObj);
//     clearClipping(alignObj);
//     AlignTraceback<TSize> trace4;
//
//     std::cout << "Score: " << globalAlignment(trace4, str, Score<int, Simple>(2, -4, -1, -2), -7, 7, BandedGotoh()) << std::endl;
//     _pumpTraceToAlign(alignObj, trace4);
//     std::cout << alignObj << std::endl;
//
//
//     clearGaps(alignObj);
//     clearClipping(alignObj);
//     AlignTraceback<TSize> trace5;
//     std::cout << "Score: " << globalAlignment(trace5, str, scoreScheme, AlignConfig<true, true, true, true>(), NeedlemanWunsch()) << std::endl;
//     _pumpTraceToAlign(alignObj, trace5);
//     std::cout << alignObj << std::endl;
//
//     clearGaps(alignObj);
//     clearClipping(alignObj);
//     AlignTraceback<TSize> trace6;
//     std::cout << "Score: " << globalAlignment(trace6, str, scoreScheme, AlignConfig<true, true, true, true>(), Gotoh()) << std::endl;
//     _pumpTraceToAlign(alignObj, trace6);
//     std::cout << alignObj << std::endl;
//     //simple trace first ... follow in general ....
//
//     //TODO (rmaerker): write test
// }

void testAlignmentTracebackImplTraceCoreFollowDiagonal()
{

    typedef Iterator<String<char> >::Type TIter;
    String<char> testMat;
    resize(testMat, 10);

    for (unsigned i = 0;i < length(testMat);++i)
    {
        testMat[i] = (char) i;
    }

    TIter it = end(testMat) - 1;

    SEQAN_ASSERT_EQ(*it, (char) 9);

    _goDiagonal(it, 5,Band<BandSwitchedOff>());

    SEQAN_ASSERT_EQ(*it, (char) 3);

    it = end(testMat) - 1;
    SEQAN_ASSERT_EQ(*it, (char) 9);
    _goDiagonal(it, 5, Band<BandSwitchedOn<> >());
    SEQAN_ASSERT_EQ(*it, (char) 4);
}

void testAlignmentTracebackImplTraceCoreFollowVertical()
{

    typedef Iterator<String<char> >::Type TIter;
    String<char> testMat;
    resize(testMat, 10);

    for (unsigned i = 0;i < length(testMat);++i)
    {
        testMat[i] = (char) i;
    }

    TIter it = end(testMat) - 1;

    SEQAN_ASSERT_EQ(*it, (char) 9);

    _goVertical(it, 5,Band<BandSwitchedOff>());

    SEQAN_ASSERT_EQ(*it, (char) 8);

    it = end(testMat) - 1;
    SEQAN_ASSERT_EQ(*it, (char) 9);
    _goVertical(it, 5, Band<BandSwitchedOn<> >());
    SEQAN_ASSERT_EQ(*it, (char) 8);
}

void testAlignmentTracebackImplTraceCoreFollowHorizontal()
{

    typedef Iterator<String<char> >::Type TIter;
    String<char> testMat;
    resize(testMat, 10);

    for (unsigned i = 0;i < length(testMat);++i)
    {
        testMat[i] = (char) i;
    }

    TIter it = end(testMat) - 2;

    SEQAN_ASSERT_EQ(*it, (char) 8);

    _goHorziontal(it, 5, Band<BandSwitchedOff>());

    SEQAN_ASSERT_EQ(*it, (char) 3);

    it = end(testMat) - 2;
    SEQAN_ASSERT_EQ(*it, (char) 8);
    _goHorziontal(it, 5, Band<BandSwitchedOn<> >());
    SEQAN_ASSERT_EQ(*it, (char) 4);
}

template <typename TGapSpec>
void testAlignmentTracebackImplTraceCoreFollowTraceLinear()
{
    String<char> seq0 = "che";
    String<char> seq1 = "cha";

    typedef AlignmentProfile<Global<>, LinearGaps, TracebackSwitchedOn > TAlignmentProfile;
    typedef typename GetTrackerSpec<TAlignmentProfile>::Type TTrackerSpec;
    typedef String<TraceBitMask::Type> TTraceMatrix;
    typedef Iterator<TTraceMatrix>::Type TIter;
    typedef Tracker<int,  TIter, Default> TTracker;



    { // Unbanded
    Band<BandSwitchedOff> band;

     TTraceMatrix traceMat;
     resize(traceMat, getTraceMatrixSize(seq0, seq1, band, True()));

     traceMat[0] = +TraceBitMask::NONE;     traceMat[4] = +TraceBitMask::HORIZONTAL; traceMat[8] = +TraceBitMask::HORIZONTAL; traceMat[12] = +TraceBitMask::HORIZONTAL;
     traceMat[1] = +TraceBitMask::VERTICAL; traceMat[5] = +TraceBitMask::VERTICAL;   traceMat[9] = +TraceBitMask::HORIZONTAL; traceMat[13] = +TraceBitMask::HORIZONTAL;
     traceMat[2] = +TraceBitMask::VERTICAL; traceMat[6] = +TraceBitMask::DIAGONAL;   traceMat[10] = +TraceBitMask::DIAGONAL;  traceMat[14] = +TraceBitMask::VERTICAL;
     traceMat[3] = +TraceBitMask::VERTICAL; traceMat[7] = +TraceBitMask::VERTICAL;   traceMat[11] = +TraceBitMask::VERTICAL;  traceMat[15] = +TraceBitMask::DIAGONAL;

    {
        TTracker tracker;
        TIter it = begin(traceMat) + 15;
        tracker(4, it);
        String<TraceSegment<size_t, size_t> > traceSegments;
        computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

        SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::DIAGONAL);
        SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 1u);
        SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 1u);
        SEQAN_ASSERT_EQ(traceSegments[0]._length, 2u);

        SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::VERTICAL);
        SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 1u);
        SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 0u);
        SEQAN_ASSERT_EQ(traceSegments[1]._length, 1u);

        SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::HORIZONTAL);
        SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 0u);
        SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 0u);
        SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);
    }


     {
         TTracker tracker;

         TIter it = begin(traceMat) + 13;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 2u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 2u);

         SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[3]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[3]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._length, 1u);
     }

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 12;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 3u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 3u);
     }

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 7;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 2u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[3]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[3]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._length, 1u);
     }

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 3;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 3u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 3u);
     }
    }

    { // Banded
     Band<BandSwitchedOn<> > band(-1,1);

     TTraceMatrix traceMat;
     resize(traceMat, getTraceMatrixSize(seq0, seq1, band, True()));

     //row 0
    //         traceMat[0] = +TraceBitMask::NONE;
     traceMat[3] = +TraceBitMask::HORIZONTAL;
     traceMat[6] = +TraceBitMask::DIAGONAL;
     traceMat[9] = (+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL);
     //row 1
     traceMat[1] = +TraceBitMask::NONE;
     traceMat[4] = +TraceBitMask::DIAGONAL;
     traceMat[7] = (+TraceBitMask::HORIZONTAL | +TraceBitMask::VERTICAL);
     traceMat[10] = (+TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL | +TraceBitMask::HORIZONTAL);
     //row 3
     traceMat[2] = +TraceBitMask::VERTICAL;
     traceMat[5] = +TraceBitMask::VERTICAL;
     traceMat[8] = (+TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL);
    //         traceMat[11] = +TraceBitMask::VERTICAL;

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 10;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

    //             for (unsigned i = 0; i < length(traceSegments);++i)
    //             {
    //                 std::cout << traceSegments[i] << std::endl;
    //             }

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[3]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[3]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._length, 1u);
     }

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 9;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 2u);

         SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);
     }

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 8;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[3]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[3]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._length, 1u);
     }
    }
}

template <typename TGapSpec>
void testAlignmentTracebackImplTraceCoreFollowTraceAffine()
{
    String<char> seq0 = "che";
    String<char> seq1 = "cha";

    typedef AlignmentProfile<Global<>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile;
    typedef typename GetTrackerSpec<TAlignmentProfile>::Type TTrackerSpec;
    typedef String<TraceBitMask::Type> TTraceMatrix;
    typedef Iterator<TTraceMatrix>::Type TIter;
    typedef Tracker<int,  TIter, Default> TTracker;



    { // Unbanded
    Band<BandSwitchedOff> band;

     TTraceMatrix traceMat;
     resize(traceMat, getTraceMatrixSize(seq0, seq1, band, True()));

     traceMat[0] = +TraceBitMask::NONE;          traceMat[4] = +TraceBitMask::HORIZONTAL_OPEN; traceMat[8] = +TraceBitMask::HORIZONTAL;      traceMat[12] = +TraceBitMask::HORIZONTAL;
     traceMat[1] = +TraceBitMask::VERTICAL_OPEN; traceMat[5] = +TraceBitMask::VERTICAL_OPEN;   traceMat[9] = +TraceBitMask::HORIZONTAL_OPEN; traceMat[13] = +TraceBitMask::HORIZONTAL;
     traceMat[2] = +TraceBitMask::VERTICAL;      traceMat[6] = +TraceBitMask::DIAGONAL;        traceMat[10] = +TraceBitMask::DIAGONAL;       traceMat[14] = +TraceBitMask::VERTICAL_OPEN;
     traceMat[3] = +TraceBitMask::VERTICAL;      traceMat[7] = +TraceBitMask::VERTICAL_OPEN;   traceMat[11] = +TraceBitMask::VERTICAL_OPEN;  traceMat[15] = +TraceBitMask::DIAGONAL;

    {
        TTracker tracker;
        TIter it = begin(traceMat) + 15;
        tracker(4, it);
        String<TraceSegment<size_t, size_t> > traceSegments;
        computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

        SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::DIAGONAL);
        SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 1u);
        SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 1u);
        SEQAN_ASSERT_EQ(traceSegments[0]._length, 2u);

        SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::VERTICAL);
        SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 1u);
        SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 0u);
        SEQAN_ASSERT_EQ(traceSegments[1]._length, 1u);

        SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::HORIZONTAL);
        SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 0u);
        SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 0u);
        SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);
    }


     {
         TTracker tracker;

         TIter it = begin(traceMat) + 13;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 2u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 2u);

         SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[3]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[3]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._length, 1u);
     }

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 12;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 3u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 3u);
     }

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 7;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 2u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[3]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[3]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._length, 1u);
     }

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 3;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 3u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 3u);
     }
    }

    { // Banded
     Band<BandSwitchedOn<> > band(-1,1);

     TTraceMatrix traceMat;
     resize(traceMat, getTraceMatrixSize(seq0, seq1, band, True()));

     //row 0
    //         traceMat[0] = +TraceBitMask::NONE;
     traceMat[3] = +TraceBitMask::HORIZONTAL_OPEN;
     traceMat[6] = +TraceBitMask::DIAGONAL;
     traceMat[9] = (+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL);
     //row 1
     traceMat[1] = +TraceBitMask::NONE;
     traceMat[4] = +TraceBitMask::DIAGONAL;
     traceMat[7] = (+TraceBitMask::HORIZONTAL_OPEN | +TraceBitMask::VERTICAL_OPEN);
     traceMat[10] = (+TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL_OPEN | +TraceBitMask::HORIZONTAL_OPEN);
     //row 3
     traceMat[2] = +TraceBitMask::VERTICAL_OPEN;
     traceMat[5] = +TraceBitMask::VERTICAL_OPEN;
     traceMat[8] = (+TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL);
    //         traceMat[11] = +TraceBitMask::VERTICAL;

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 10;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

    //             for (unsigned i = 0; i < length(traceSegments);++i)
    //             {
    //                 std::cout << traceSegments[i] << std::endl;
    //             }

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[3]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[3]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._length, 1u);
     }

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 9;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 2u);

         SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[3]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[3]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._length, 1u);
     }

     {
         TTracker tracker;

         TIter it = begin(traceMat) + 8;
         tracker(4, it);
         String<TraceSegment<size_t, size_t> > traceSegments;
         computeTraceback(traceSegments, tracker, traceMat, seq0, seq1, band, TAlignmentProfile());

         SEQAN_ASSERT_EQ(traceSegments[0]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[0]._horizontalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[0]._verticalBeginPos, 3u);
         SEQAN_ASSERT_EQ(traceSegments[0]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[1]._traceValue, +TraceBitMask::VERTICAL);
         SEQAN_ASSERT_EQ(traceSegments[1]._horizontalBeginPos, 2u);
         SEQAN_ASSERT_EQ(traceSegments[1]._verticalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[1]._length, 2u);

         SEQAN_ASSERT_EQ(traceSegments[2]._traceValue, +TraceBitMask::DIAGONAL);
         SEQAN_ASSERT_EQ(traceSegments[2]._horizontalBeginPos, 1u);
         SEQAN_ASSERT_EQ(traceSegments[2]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[2]._length, 1u);

         SEQAN_ASSERT_EQ(traceSegments[3]._traceValue, +TraceBitMask::HORIZONTAL);
         SEQAN_ASSERT_EQ(traceSegments[3]._horizontalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._verticalBeginPos, 0u);
         SEQAN_ASSERT_EQ(traceSegments[3]._length, 1u);
     }
    }
}

SEQAN_DEFINE_TEST(test_alignment_traceback_impl_get_matrix_size)
{
    testAlignmentTracebackImplTraceCoreGetTraceMatrixSize();
}


SEQAN_DEFINE_TEST(test_alignment_traceback_impl_follow_diagonal)
{
    testAlignmentTracebackImplTraceCoreFollowDiagonal();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_impl_follow_vertical)
{
    testAlignmentTracebackImplTraceCoreFollowVertical();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_impl_follow_horizontal)
{
    testAlignmentTracebackImplTraceCoreFollowHorizontal();
}

SEQAN_DEFINE_TEST(test_alignment_traceback_impl_follow_trace)
{
     testAlignmentTracebackImplTraceCoreFollowTraceLinear<LinearGaps>();
     testAlignmentTracebackImplTraceCoreFollowTraceAffine<AffineGaps>();
}

#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_TRACEBACK_IMPL_H_
