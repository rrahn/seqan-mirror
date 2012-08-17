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

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_TRACKER_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_TRACKER_H_

#include <seqan/basic.h>

#include <seqan/align/alignment_base.h>
#include <seqan/align/alignment_dp_band.h>
#include <seqan/align/alignment_dp_manager.h>
#include <seqan/align/alignment_dp_value.h>
#include <seqan/align/alignment_dp_tracker.h>

using namespace seqan;


void
testAlignmentDPTrackerGridCoordinate()
{
    GridCoordinate coordinate1;

    SEQAN_ASSERT_EQ(coordinate1.x, 0u);
    SEQAN_ASSERT_EQ(coordinate1.y, 0u);

    size_t _x = 10;
    size_t _y = 112;
    GridCoordinate  coordinate2(_x,_y);

    SEQAN_ASSERT_EQ(coordinate2.x, 10u);
    SEQAN_ASSERT_EQ(coordinate2.y, 112u);

    coordinate2.x = 13u;
    coordinate2.y = 13u;

    SEQAN_ASSERT_EQ(coordinate2.x, 13u);
    SEQAN_ASSERT_EQ(coordinate2.y, 13u);

    coordinate2.x = -12;
    coordinate2.y = 45;

    SEQAN_ASSERT_GEQ(coordinate2.x, 0u);
    SEQAN_ASSERT_EQ(coordinate2.y, 45u);
}

void
testAlignmentDPTrackerGetTrackerSpec()
{
    { // no traceback
        typedef AlignmentProfile<Global<>, DPValue<int, LinearGaps>, TracebackSwitchedOff > TProfile;
        bool result = IsSameType<GetTrackerSpec<TProfile>::Type, ScoreOnly>::VALUE;
        SEQAN_ASSERT_EQ(result, true);
    }

    { // default traceback
        typedef AlignmentProfile<Global<>, DPValue<int, LinearGaps>, TracebackSwitchedOn > TProfile;
        bool result = IsSameType<GetTrackerSpec<TProfile>::Type, Default>::VALUE;
        SEQAN_ASSERT_EQ(result, true);
    }

    { // default traceback
        typedef AlignmentProfile<Global<FreeEndGaps<> >, DPValue<int, LinearGaps>, TracebackSwitchedOn > TProfile;
        bool result = IsSameType<GetTrackerSpec<TProfile>::Type, Default>::VALUE;
        SEQAN_ASSERT_EQ(result, true);
    }

    { // default traceback
        typedef AlignmentProfile<Local<>, DPValue<int, LinearGaps>, TracebackSwitchedOn > TProfile;
        bool result = IsSameType<GetTrackerSpec<TProfile>::Type, Default>::VALUE;
        SEQAN_ASSERT_EQ(result, true);
    }

    { // banded chain traceback
        typedef AlignmentProfile<Global<BandedChain>, DPValue<int, LinearGaps>, TracebackSwitchedOn > TProfile;
        bool result = IsSameType<GetTrackerSpec<TProfile>::Type, BandedChain>::VALUE;
        SEQAN_ASSERT_EQ(result, true);
    }

    { // row max traceback
        typedef AlignmentProfile<Global<RowMax>, DPValue<int, LinearGaps>, TracebackSwitchedOn > TProfile;
        bool result = IsSameType<GetTrackerSpec<TProfile>::Type, RowMax>::VALUE;
        SEQAN_ASSERT_EQ(result, true);
    }

    { // waterman eggert traceback
        typedef AlignmentProfile<Local<WatermanEggert>, DPValue<int, LinearGaps>, TracebackSwitchedOn > TProfile;
        bool result = IsSameType<GetTrackerSpec<TProfile>::Type, WatermanEggert>::VALUE;
        SEQAN_ASSERT_EQ(result, true);
    }
}

void
testAlignmentDPTrackerScoreOnly()
{
    typedef AlignmentProfile<Global<>, DPValue<int, LinearGaps>, TracebackSwitchedOff > TProfile;

    typedef Tracker<int, CharIterator, GetTrackerSpec<TProfile>::Type > TTracker;

    CharString seq = "This is a test";
    CharIterator it = begin(seq) + 10;



    TTracker tracker;
    SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);

    tracker(10, it);

    SEQAN_ASSERT_EQ(tracker._maxScore, 10);
    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);

    tracker(20, ++it);
    SEQAN_ASSERT_EQ(tracker._maxScore, 20);
    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);

    tracker(-12, ++it);
    SEQAN_ASSERT_EQ(tracker._maxScore, 20);
    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);
}

void
testAlignmentDPTrackerDefault()
{
    typedef AlignmentProfile<Global<>, DPValue<int, LinearGaps>, TracebackSwitchedOn > TProfile;

    typedef Tracker<int, CharIterator, GetTrackerSpec<TProfile>::Type > TTracker;
                    //01234567890
    CharString seq = "This is a test";
    CharIterator it = begin(seq) + 10;

    TTracker tracker;
    SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 1u);

    tracker(10, it);

    SEQAN_ASSERT_EQ(tracker._maxScore, 10);
    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 1u);
    SEQAN_ASSERT_EQ(value(tracker._tracePoints[0]), 't');

    tracker(20, ++it);
    SEQAN_ASSERT_EQ(tracker._maxScore, 20);
    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 1u);
    SEQAN_ASSERT_EQ(value(tracker._tracePoints[0]), 'e');

    tracker(-12, ++it);
    SEQAN_ASSERT_EQ(tracker._maxScore, 20);
    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 1u);
    SEQAN_ASSERT_EQ(value(tracker._tracePoints[0]), 'e');
}

void
testAlignmentDPTrackerRowMax()
{
    typedef AlignmentProfile<Global<RowMax>, DPValue<int, LinearGaps>, TracebackSwitchedOn > TProfile;

    typedef Tracker<int, CharIterator, GetTrackerSpec<TProfile>::Type > TTracker;
                    //01234567890
    CharString seq = "This is a test";
    CharIterator it = begin(seq) + 10;

    SEQAN_FAIL("TODO: Implement this test");
//  TTracker tracker;
//  SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
//  SEQAN_ASSERT_EQ(length(tracker._tracePoints), 1u);
//
//  tracker(10, it);
//
//  SEQAN_ASSERT_EQ(tracker._maxScore, 10);
//  SEQAN_ASSERT_EQ(length(tracker._tracePoints), 1u);
//  SEQAN_ASSERT_EQ(value(tracker._tracePoints[0]), 't');
//
//  tracker(20, ++it);
//  SEQAN_ASSERT_EQ(tracker._maxScore, 20);
//  SEQAN_ASSERT_EQ(length(tracker._tracePoints), 1u);
//  SEQAN_ASSERT_EQ(value(tracker._tracePoints[0]), 'e');
//
//  tracker(-12, ++it);
//  SEQAN_ASSERT_EQ(tracker._maxScore, 20);
//  SEQAN_ASSERT_EQ(length(tracker._tracePoints), 1u);
//  SEQAN_ASSERT_EQ(value(tracker._tracePoints[0]), 'e');
}

void
testAlignmentDPTrackerBandedChainInit()
{
//    typedef AlignmentProfile<Global<BandedChain>, DPValue<int, AffineGaps>, Traceback<On> > TProfile;
//    typedef DPValue<int, AffineGaps> TDPValue;
//    typedef Tracker<int, CharIterator, GetTrackerSpec<TProfile>::Type > TTracker;
//
//                    //01234567890123456
//    CharString mat = "actgactgacgatcgac";
//    // need a different test ...
//    CharIterator it = begin(mat);
//
//    TTracker tracker;
//    SEQAN_ASSERT_EQ(tracker.x, 0u);
//    SEQAN_ASSERT_EQ(tracker.y, 0u);
//    SEQAN_ASSERT_EQ(tracker._posDim0, 0u);
//    SEQAN_ASSERT_EQ(tracker._posDim1, 0u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim0Curr), 0u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim1Curr), 0u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim0Next), 0u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim1Next), 0u);
//    SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
//    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);
//
//    _init(tracker,4u,4u,2,2,2u,2u);
//    SEQAN_ASSERT_EQ(tracker.x, 0u);
//    SEQAN_ASSERT_EQ(tracker.y, 0u);
//    SEQAN_ASSERT_EQ(tracker._posDim0, 2u);
//    SEQAN_ASSERT_EQ(tracker._posDim1, 2u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim0Curr), 4u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim1Curr), 4u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim0Next), 2u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim1Next), 2u);
//    SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
//    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);
//
//    SEQAN_ASSERT_EQ(tracker._initDim0Curr[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim0Curr[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim0Curr[2], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim0Curr[3], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Curr[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Curr[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Curr[2], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Curr[3], TDPValue(0,+TraceBitMask::FORBIDDEN));
//
//    SEQAN_ASSERT_EQ(tracker._initDim0Next[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim0Next[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Next[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Next[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//
//    tracker._initDim0Curr[2] = TDPValue(20,+TraceBitMask::DIAGONAL);
//    tracker._initDim1Curr[0] = TDPValue(20,+TraceBitMask::DIAGONAL);
//    tracker._initDim0Next[0] = TDPValue(20,+TraceBitMask::DIAGONAL);
//    tracker._initDim1Next[1] = TDPValue(20,+TraceBitMask::VERTICAL);
//    tracker.x = 12;
//    tracker.y = 12;
//
//    _init(tracker, 2u, 2u, 3,3, 4u,4u);
//
//    SEQAN_ASSERT_EQ(tracker.x, 0u);
//    SEQAN_ASSERT_EQ(tracker.y, 0u);
//    SEQAN_ASSERT_EQ(tracker._posDim0, 3u);
//    SEQAN_ASSERT_EQ(tracker._posDim1, 3u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim0Curr), 4u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim1Curr), 4u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim0Next), 4u);
//    SEQAN_ASSERT_EQ(length(tracker._initDim1Next), 4u);
//    SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
//    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);
//
//    SEQAN_ASSERT_EQ(tracker._initDim0Curr[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim0Curr[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim0Curr[2], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim0Curr[3], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Curr[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Curr[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Curr[2], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Curr[3], TDPValue(0,+TraceBitMask::FORBIDDEN));
//
//    SEQAN_ASSERT_EQ(tracker._initDim0Next[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim0Next[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim0Next[2], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim0Next[3], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Next[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Next[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Next[2], TDPValue(0,+TraceBitMask::FORBIDDEN));
//    SEQAN_ASSERT_EQ(tracker._initDim1Next[3], TDPValue(0,+TraceBitMask::FORBIDDEN));
}

void
testAlignmentDPTrackerGoBeginNextColumn()
{
    typedef AlignmentProfile<Global<BandedChain>, DPValue<int, AffineGaps>, TracebackSwitchedOn > TProfile;
    typedef DPValue<int, AffineGaps> TDPValue;
    typedef Tracker<int, CharIterator, GetTrackerSpec<TProfile>::Type > TTracker;

                    //01234567890123456
    CharString mat = "actgactgacgatcgac";
    // need a different test ...
    CharIterator it = begin(mat);

    TTracker tracker;
    SEQAN_ASSERT_EQ(tracker.x, 0u);
    SEQAN_ASSERT_EQ(tracker.y, 0u);
    SEQAN_ASSERT_EQ(tracker._posDim0, 0u);
    SEQAN_ASSERT_EQ(tracker._posDim1, 0u);
    SEQAN_ASSERT_EQ(length(tracker._initDim0Curr), 0u);
    SEQAN_ASSERT_EQ(length(tracker._initDim1Curr), 0u);
    SEQAN_ASSERT_EQ(length(tracker._initDim0Next), 0u);
    SEQAN_ASSERT_EQ(length(tracker._initDim1Next), 0u);
    SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);

    goBeginNextColumn(tracker, -1, Band<BandSwitchedOn<> >(2,2));
    SEQAN_ASSERT_EQ(tracker.x, 1u);
    SEQAN_ASSERT_EQ(tracker.y, 0u);
    goBeginNextColumn(tracker, 0, Band<BandSwitchedOn<> >(2,2));
    SEQAN_ASSERT_EQ(tracker.x, 2u);
    SEQAN_ASSERT_EQ(tracker.y, 1u);
    goBeginNextColumn(tracker, 1, Band<BandSwitchedOn<> >(2,2));
    SEQAN_ASSERT_EQ(tracker.x, 3u);
    SEQAN_ASSERT_EQ(tracker.y, 2u);
    goBeginNextColumn(tracker, 2, Band<BandSwitchedOn<> >(2,2));
    SEQAN_ASSERT_EQ(tracker.x, 4u);
    SEQAN_ASSERT_EQ(tracker.y, 3u);

    goBeginNextColumn(tracker, 0, Band<BandSwitchedOff>());
    SEQAN_ASSERT_EQ(tracker.x, 5u);
    SEQAN_ASSERT_EQ(tracker.y, 0u);
    goBeginNextColumn(tracker, 1, Band<BandSwitchedOff>());
    SEQAN_ASSERT_EQ(tracker.x, 6u);
    SEQAN_ASSERT_EQ(tracker.y, 0u);
    goBeginNextColumn(tracker, 2, Band<BandSwitchedOff>());
    SEQAN_ASSERT_EQ(tracker.x, 7u);
    SEQAN_ASSERT_EQ(tracker.y, 0u);
    goBeginNextColumn(tracker, 3, Band<BandSwitchedOff>());
    SEQAN_ASSERT_EQ(tracker.x, 8u);
    SEQAN_ASSERT_EQ(tracker.y, 0u);
}

void
testAlignmentDPTrackerBandedChain()
{
//    typedef AlignmentProfile<Global<BandedChain>, DPValue<int, AffineGaps>, Traceback<On> > TProfile;
//    typedef DPValue<int, AffineGaps> TDPValue;
//    typedef Tracker<int, CharIterator, GetTrackerSpec<TProfile>::Type > TTracker;
//
//                    //0123456789012345
//    CharString mat = "actgactgacgatcga";
//    // need a different test ...
//    CharIterator it = begin(mat);
//
//    { // default c'tor;
//        TTracker tracker;
//        SEQAN_ASSERT_EQ(tracker.x, 0u);
//        SEQAN_ASSERT_EQ(tracker.y, 0u);
//        SEQAN_ASSERT_EQ(tracker._posDim0, 0u);
//        SEQAN_ASSERT_EQ(tracker._posDim1, 0u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim0Curr), 0u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim1Curr), 0u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim0Next), 0u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim1Next), 0u);
//        SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
//        SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);
//
//    }
//
//    { // c'tor and modification
//        TTracker tracker(2u,2u,2,2,2u,2u);
//        SEQAN_ASSERT_EQ(tracker.x, 0u);
//        SEQAN_ASSERT_EQ(tracker.y, 0u);
//        SEQAN_ASSERT_EQ(tracker._posDim0, 2u);
//        SEQAN_ASSERT_EQ(tracker._posDim1, 2u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim0Curr), 2u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim1Curr), 2u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim0Next), 2u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim1Next), 2u);
//        SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
//        SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);
//
//        SEQAN_ASSERT_EQ(tracker._initDim0Curr[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim0Curr[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim1Curr[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim1Curr[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//
//        SEQAN_ASSERT_EQ(tracker._initDim0Next[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim0Next[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim1Next[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim1Next[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//
//        //first column
//        tracker(0, it); tracker(TDPValue(0, +TraceBitMask::NONE));
//        SEQAN_ASSERT_EQ(tracker.x, 0u);
//        SEQAN_ASSERT_EQ(tracker.y, 1u);
//        tracker(1, ++it); tracker(TDPValue(1, +TraceBitMask::VERTICAL));
//        SEQAN_ASSERT_EQ(tracker.x, 0u);
//        SEQAN_ASSERT_EQ(tracker.y, 2u);
//        tracker(-12, ++it); tracker(TDPValue(-12, +TraceBitMask::VERTICAL));
//        SEQAN_ASSERT_EQ(tracker.x, 0u);
//        SEQAN_ASSERT_EQ(tracker.y, 3u);
//        tracker(14, ++it); tracker(TDPValue(14, +TraceBitMask::VERTICAL));
//        SEQAN_ASSERT_EQ(tracker.x, 0u);
//        SEQAN_ASSERT_EQ(tracker.y, 4u);
//
//        goBeginNextColumn(tracker, 0, Band<BandSwitchedOff>());
//
//        //second column
//        SEQAN_ASSERT_EQ(tracker.x, 1u);
//        SEQAN_ASSERT_EQ(tracker.y, 0u);
//        tracker(14, ++it); tracker(TDPValue(14, +TraceBitMask::HORIZONTAL));
//        SEQAN_ASSERT_EQ(tracker.x, 1u);
//        SEQAN_ASSERT_EQ(tracker.y, 1u);
//        tracker(1, ++it); tracker(TDPValue(1, +TraceBitMask::DIAGONAL));
//        SEQAN_ASSERT_EQ(tracker.x, 1u);
//        SEQAN_ASSERT_EQ(tracker.y, 2u);
//        tracker(0, ++it);tracker(TDPValue(0, +TraceBitMask::VERTICAL));
//        SEQAN_ASSERT_EQ(tracker.x, 1u);
//        SEQAN_ASSERT_EQ(tracker.y, 3u);
//        tracker(1, ++it); tracker(TDPValue(1, +TraceBitMask::VERTICAL));
//        SEQAN_ASSERT_EQ(tracker.x, 1u);
//        SEQAN_ASSERT_EQ(tracker.y, 4u);
//
//        goBeginNextColumn(tracker, 0, Band<BandSwitchedOff>());
//
//        //third column
//        SEQAN_ASSERT_EQ(tracker.x, 2u);
//        SEQAN_ASSERT_EQ(tracker.y, 0u);
//        tracker(0, ++it); tracker(TDPValue(0, +TraceBitMask::HORIZONTAL));
//        SEQAN_ASSERT_EQ(tracker.x, 2u);
//        SEQAN_ASSERT_EQ(tracker.y, 1u);
//        tracker(1, ++it); tracker(TDPValue(1, +TraceBitMask::HORIZONTAL));
//        SEQAN_ASSERT_EQ(tracker.x, 2u);
//        SEQAN_ASSERT_EQ(tracker.y, 2u);
//        tracker(20, ++it); tracker(TDPValue(20, +TraceBitMask::DIAGONAL));
//        SEQAN_ASSERT_EQ(tracker.x, 2u);
//        SEQAN_ASSERT_EQ(tracker.y, 3u);
//        tracker(1, ++it); tracker(TDPValue(1, +TraceBitMask::VERTICAL));
//        SEQAN_ASSERT_EQ(tracker.x, 2u);
//        SEQAN_ASSERT_EQ(tracker.y, 4u);
//
//        goBeginNextColumn(tracker, 0, Band<BandSwitchedOff>());
//
//        //fourth column
//        SEQAN_ASSERT_EQ(tracker.x, 3u);
//        SEQAN_ASSERT_EQ(tracker.y, 0u);
//        tracker(0, ++it); tracker(TDPValue(0, +TraceBitMask::HORIZONTAL));
//        SEQAN_ASSERT_EQ(tracker.x, 3u);
//        SEQAN_ASSERT_EQ(tracker.y, 1u);
//        tracker(1, ++it); tracker(TDPValue(1, +TraceBitMask::HORIZONTAL));
//        SEQAN_ASSERT_EQ(tracker.x, 3u);
//        SEQAN_ASSERT_EQ(tracker.y, 2u);
//        tracker(0, ++it); tracker(TDPValue(0, +TraceBitMask::HORIZONTAL));
//        SEQAN_ASSERT_EQ(tracker.x, 3u);
//        SEQAN_ASSERT_EQ(tracker.y, 3u);
//        tracker(20, ++it); tracker(TDPValue(20, +TraceBitMask::DIAGONAL));
//        SEQAN_ASSERT_EQ(tracker.x, 3u);
//        SEQAN_ASSERT_EQ(tracker.y, 4u);
//
//        SEQAN_ASSERT_EQ(tracker._initDim0Next[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim0Next[1], TDPValue(0,+TraceBitMask::HORIZONTAL));
//
//        SEQAN_ASSERT_EQ(tracker._initDim1Next[0], TDPValue(20,+TraceBitMask::DIAGONAL));
//        SEQAN_ASSERT_EQ(tracker._initDim1Next[1], TDPValue(1,+TraceBitMask::VERTICAL));
//        SEQAN_ASSERT_EQ(tracker._maxScore, 20);
//        SEQAN_ASSERT_EQ(length(tracker._tracePoints), 2u);
//        SEQAN_ASSERT_EQ(value(tracker._tracePoints[0]), 'g');
//        SEQAN_ASSERT_EQ(value(tracker._tracePoints[1]), 'a');
//    }
}

void
testAlignmentDPTrackerBandedChainBanded()
{
//    typedef AlignmentProfile<Global<BandedChain>, DPValue<int, AffineGaps>, Traceback<On> > TProfile;
//    typedef DPValue<int, AffineGaps> TDPValue;
//    typedef Tracker<int, CharIterator, GetTrackerSpec<TProfile>::Type > TTracker;
//
//                    //0123456789012345
//    CharString mat = "actgactgacgatcga";
//    // need a different test ...
//    CharIterator it = begin(mat);
//
//    { // default c'tor;
//        TTracker tracker;
//        SEQAN_ASSERT_EQ(tracker.x, 0u);
//        SEQAN_ASSERT_EQ(tracker.y, 0u);
//        SEQAN_ASSERT_EQ(tracker._posDim0, 0u);
//        SEQAN_ASSERT_EQ(tracker._posDim1, 0u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim0Curr), 0u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim1Curr), 0u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim0Next), 0u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim1Next), 0u);
//        SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
//        SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);
//
//    }
//
//    { // c'tor and modification
//        TTracker tracker(2u,2u,2,2,2u,2u);
//        SEQAN_ASSERT_EQ(tracker.x, 0u);
//        SEQAN_ASSERT_EQ(tracker.y, 0u);
//        SEQAN_ASSERT_EQ(tracker._posDim0, 2u);
//        SEQAN_ASSERT_EQ(tracker._posDim1, 2u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim0Curr), 2u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim1Curr), 2u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim0Next), 2u);
//        SEQAN_ASSERT_EQ(length(tracker._initDim1Next), 2u);
//        SEQAN_ASSERT_EQ(tracker._maxScore, MinValue<int>::VALUE);
//        SEQAN_ASSERT_EQ(length(tracker._tracePoints), 0u);
//
//        SEQAN_ASSERT_EQ(tracker._initDim0Curr[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim0Curr[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim1Curr[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim1Curr[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//
//        SEQAN_ASSERT_EQ(tracker._initDim0Next[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim0Next[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim1Next[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim1Next[1], TDPValue(0,+TraceBitMask::FORBIDDEN));
//
//        //first column
//        tracker(0, it); tracker(TDPValue(0, +TraceBitMask::NONE));
//        SEQAN_ASSERT_EQ(tracker.x, 0u);
//        SEQAN_ASSERT_EQ(tracker.y, 1u);
//        tracker(1, ++it); tracker(TDPValue(1, +TraceBitMask::VERTICAL));
//        SEQAN_ASSERT_EQ(tracker.x, 0u);
//        SEQAN_ASSERT_EQ(tracker.y, 2u);
//
//        goBeginNextColumn(tracker, -1, Band<int, On>(5,6));
//        it += 2;
//
//        //second column
//        SEQAN_ASSERT_EQ(tracker.x, 1u);
//        SEQAN_ASSERT_EQ(tracker.y, 0u);
//        tracker(14, ++it); tracker(TDPValue(14, +TraceBitMask::HORIZONTAL));
//        SEQAN_ASSERT_EQ(tracker.x, 1u);
//        SEQAN_ASSERT_EQ(tracker.y, 1u);
//        tracker(1, ++it); tracker(TDPValue(1, +TraceBitMask::DIAGONAL));
//        SEQAN_ASSERT_EQ(tracker.x, 1u);
//        SEQAN_ASSERT_EQ(tracker.y, 2u);
//        tracker(0, ++it);tracker(TDPValue(0, +TraceBitMask::VERTICAL));
//        SEQAN_ASSERT_EQ(tracker.x, 1u);
//        SEQAN_ASSERT_EQ(tracker.y, 3u);
//
//        goBeginNextColumn(tracker, 0, Band<int, On>(5,6));
//        it += 2;
//
//        //third column
//        SEQAN_ASSERT_EQ(tracker.x, 2u);
//        SEQAN_ASSERT_EQ(tracker.y, 1u);
//        tracker(1, ++it); tracker(TDPValue(1, +TraceBitMask::HORIZONTAL));
//        SEQAN_ASSERT_EQ(tracker.x, 2u);
//        SEQAN_ASSERT_EQ(tracker.y, 2u);
//        tracker(20, ++it); tracker(TDPValue(20, +TraceBitMask::DIAGONAL));
//        SEQAN_ASSERT_EQ(tracker.x, 2u);
//        SEQAN_ASSERT_EQ(tracker.y, 3u);
//        tracker(1, ++it); tracker(TDPValue(1, +TraceBitMask::VERTICAL));
//        SEQAN_ASSERT_EQ(tracker.x, 2u);
//        SEQAN_ASSERT_EQ(tracker.y, 4u);
//
//        goBeginNextColumn(tracker, 1, Band<int, On>(5,6));
//        it += 2;
//
//        //fourth column
//        SEQAN_ASSERT_EQ(tracker.x, 3u);
//        SEQAN_ASSERT_EQ(tracker.y, 2u);
//        tracker(0, ++it); tracker(TDPValue(0, +TraceBitMask::HORIZONTAL));
//        SEQAN_ASSERT_EQ(tracker.x, 3u);
//        SEQAN_ASSERT_EQ(tracker.y, 3u);
//        tracker(20, ++it); tracker(TDPValue(20, +TraceBitMask::DIAGONAL));
//        SEQAN_ASSERT_EQ(tracker.x, 3u);
//        SEQAN_ASSERT_EQ(tracker.y, 4u);
//
//        SEQAN_ASSERT_EQ(tracker._initDim0Next[0], TDPValue(0,+TraceBitMask::FORBIDDEN));
//        SEQAN_ASSERT_EQ(tracker._initDim0Next[1], TDPValue(0,+TraceBitMask::HORIZONTAL));
//
//        SEQAN_ASSERT_EQ(tracker._initDim1Next[0], TDPValue(20,+TraceBitMask::DIAGONAL));
//        SEQAN_ASSERT_EQ(tracker._initDim1Next[1], TDPValue(1,+TraceBitMask::VERTICAL));
//
//        SEQAN_ASSERT_EQ(tracker._maxScore, 20);
//        SEQAN_ASSERT_EQ(length(tracker._tracePoints), 2u);
//        SEQAN_ASSERT_EQ(value(tracker._tracePoints[0]), 'g');
//        SEQAN_ASSERT_EQ(value(tracker._tracePoints[1]), 'a');
//    }
}

SEQAN_DEFINE_TEST(test_alignment_dp_tracker_grid_coordinate)
{
    testAlignmentDPTrackerGridCoordinate();
}

SEQAN_DEFINE_TEST(test_alignment_dp_tracker_get_tracker_spec)
{
    testAlignmentDPTrackerGetTrackerSpec();
}

SEQAN_DEFINE_TEST(test_alignment_dp_tracker_score_only)
{
    testAlignmentDPTrackerScoreOnly();
}

SEQAN_DEFINE_TEST(test_alignment_dp_tracker_default)
{
    testAlignmentDPTrackerDefault();
}

SEQAN_DEFINE_TEST(test_alignment_dp_tracker_row_max)
{
    SEQAN_FAIL("TODO: Implement test");
}

SEQAN_DEFINE_TEST(test_alignment_dp_tracker_waterman_eggert)
{
    SEQAN_FAIL("TODO: Implement test");
}

SEQAN_DEFINE_TEST(test_alignment_dp_tracker_banded_chain_init)
{
    testAlignmentDPTrackerBandedChainInit();
}

SEQAN_DEFINE_TEST(test_alignment_dp_tracker_go_begin_next_column)
{
    testAlignmentDPTrackerGoBeginNextColumn();
}

SEQAN_DEFINE_TEST(test_alignment_dp_tracker_banded_chain)
{
    testAlignmentDPTrackerBandedChain();
    testAlignmentDPTrackerBandedChainBanded();
}

#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_TRACKER_H_
