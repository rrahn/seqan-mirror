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

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_MANAGER_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_MANAGER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/align/alignment_base.h>
#include <seqan/align/alignment_dp_value.h>
#include <seqan/align/alignment_dp_band.h>
#include <seqan/align/alignment_dp_manager.h>

using namespace seqan;


void testAlignmentManagerGetDpDirection()
{
    typedef Cell<DPDirectionAll> TCell1;
    typedef IsSameType<DPDirectionAll, GetDpDirection<TCell1>::Type>::Type TResult1;

    SEQAN_ASSERT_EQ(+TResult1::VALUE, true);

    typedef Cell<DPDirectionUpperBand> TCell2;
    typedef IsSameType<DPDirectionUpperBand, GetDpDirection<TCell2>::Type>::Type TResult2;

    SEQAN_ASSERT_EQ(+TResult2::VALUE, true);
}

void testAlignmentManagerIsToTrack()
{
    typedef Cell<DPDirectionAll> TCell1;
    SEQAN_ASSERT_EQ(+IsToTrack<TCell1>::VALUE, false);

    typedef Cell<DPDirectionUpperBand, True> TCell2;
    SEQAN_ASSERT_EQ(+IsToTrack<TCell2>::VALUE, true);

    typedef Cell<DPDirectionAll, False> TCell3;
    SEQAN_ASSERT_EQ(+IsToTrack<TCell3>::VALUE, false);
}


void testAlignmentManagerSetUpColumnManager()
{

    typedef Band<BandSwitchedOff> TBand1;
    typedef Band<BandSwitchedOn<> > TBand2;

    typedef AlignmentProfile<Global<FreeEndGaps<> >, LinearGaps, TracebackSwitchedOff > TAligmentProfile1;
    typedef AlignmentProfile<Global<FreeEndGaps<True, True, True, True> >, LinearGaps, TracebackSwitchedOn > TAligmentProfile2;
    typedef AlignmentProfile<Local<>, LinearGaps, TracebackSwitchedOff > TAligmentProfile3;
    typedef AlignmentProfile<Global<SplitBreakpoint>, LinearGaps, TracebackSwitchedOn > TAligmentProfile4;

    { // standard global and unbanded
        typedef SetUpColumnManager<TAligmentProfile1, DPNoBandInitPhase> TCM1;

        typedef IsSameType<TCM1::TFirstCell_, Cell<DPDirectionZero, False> >::Type TResult1a;
        typedef IsSameType<TCM1::TInnerCell_, Cell<DPDirectionVertical, False> >::Type TResult1b;
        typedef IsSameType<TCM1::TLastCell_, Cell<DPDirectionVertical, False> >::Type TResult1c;

        SEQAN_ASSERT_EQ(+TResult1a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1c::VALUE, true);


        // does not exists anymore
//        typedef SetUpColumnManager<TAligmentProfile1, TBand1, BandOpenColumn> TCM2;
//
//        typedef IsSameType<TCM2::TFirstCell_, Cell<DPDirectionHorizontal, False> >::Type TResult2a;
//        typedef IsSameType<TCM2::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult2b;
//        typedef IsSameType<TCM2::TLastCell_, Cell<DPDirectionLowerBand, False> >::Type TResult2c;
//
//        SEQAN_ASSERT_EQ(+TResult2a::VALUE, true);
//        SEQAN_ASSERT_EQ(+TResult2b::VALUE, true);
//        SEQAN_ASSERT_EQ(+TResult2c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile1, DPNoBandPhase> TCM3;

        typedef IsSameType<TCM3::TFirstCell_, Cell<DPDirectionHorizontal, False> >::Type TResult3a;
        typedef IsSameType<TCM3::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult3b;
        typedef IsSameType<TCM3::TLastCell_, Cell<DPDirectionAll, False> >::Type TResult3c;

        SEQAN_ASSERT_EQ(+TResult3a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3c::VALUE, true);

//        typedef SetUpColumnManager<TAligmentProfile1, TBand1, BandCloseColumn> TCM4;
//
//        typedef IsSameType<TCM4::TFirstCell_, Cell<DPDirectionUpperBand, False> >::Type TResult4a;
//        typedef IsSameType<TCM4::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult4b;
//        typedef IsSameType<TCM4::TLastCell_, Cell<DPDirectionAll, False> >::Type TResult4c;
//
//        SEQAN_ASSERT_EQ(+TResult4a::VALUE, true);
//        SEQAN_ASSERT_EQ(+TResult4b::VALUE, true);
//        SEQAN_ASSERT_EQ(+TResult4c::VALUE, true);
    }

    {   // standard global and banded
        typedef SetUpColumnManager<TAligmentProfile1, DPBandInitPhase> TCM1;

        typedef IsSameType<TCM1::TFirstCell_, Cell<DPDirectionZero, False> >::Type TResult1a;
        typedef IsSameType<TCM1::TInnerCell_, Cell<DPDirectionVertical, False> >::Type TResult1b;
        typedef IsSameType<TCM1::TLastCell_, Cell<DPDirectionVertical, False> >::Type TResult1c;

        SEQAN_ASSERT_EQ(+TResult1a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile1, DPBandFirstPhase> TCM2;

        typedef IsSameType<TCM2::TFirstCell_, Cell<DPDirectionHorizontal, False> >::Type TResult2a;
        typedef IsSameType<TCM2::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult2b;
        typedef IsSameType<TCM2::TLastCell_, Cell<DPDirectionLowerBand, False> >::Type TResult2c;

        SEQAN_ASSERT_EQ(+TResult2a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult2b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult2c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile1, DPBandMiddlePhase> TCM3;

        typedef IsSameType<TCM3::TFirstCell_, Cell<DPDirectionUpperBand, False> >::Type TResult3a;
        typedef IsSameType<TCM3::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult3b;
        typedef IsSameType<TCM3::TLastCell_, Cell<DPDirectionLowerBand, False> >::Type TResult3c;

        SEQAN_ASSERT_EQ(+TResult3a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile1, DPBandLastPhase> TCM4;

        typedef IsSameType<TCM4::TFirstCell_, Cell<DPDirectionUpperBand, False> >::Type TResult4a;
        typedef IsSameType<TCM4::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult4b;
        typedef IsSameType<TCM4::TLastCell_, Cell<DPDirectionAll, False> >::Type TResult4c;

        SEQAN_ASSERT_EQ(+TResult4a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult4b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult4c::VALUE, true);
    }

    {   // overlap and unbanded
        typedef SetUpColumnManager<TAligmentProfile2, DPNoBandInitPhase> TCM1;

        typedef IsSameType<TCM1::TFirstCell_, Cell<DPDirectionZero, False> >::Type TResult1a;
        typedef IsSameType<TCM1::TInnerCell_, Cell<DPDirectionZero, False> >::Type TResult1b;
        typedef IsSameType<TCM1::TLastCell_, Cell<DPDirectionZero, True> >::Type TResult1c;

        SEQAN_ASSERT_EQ(+TResult1a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1c::VALUE, true);


        typedef SetUpColumnManager<TAligmentProfile2, DPNoBandPhase> TCM3;

        typedef IsSameType<TCM3::TFirstCell_, Cell<DPDirectionZero, False> >::Type TResult3a;
        typedef IsSameType<TCM3::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult3b;
        typedef IsSameType<TCM3::TLastCell_, Cell<DPDirectionAll, True> >::Type TResult3c;

        SEQAN_ASSERT_EQ(+TResult3a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3c::VALUE, true);
    }

    {   // overlap and banded
        typedef SetUpColumnManager<TAligmentProfile2, DPBandInitPhase> TCM1;

        typedef IsSameType<TCM1::TFirstCell_, Cell<DPDirectionZero, False> >::Type TResult1a;
        typedef IsSameType<TCM1::TInnerCell_, Cell<DPDirectionZero, False> >::Type TResult1b;
        typedef IsSameType<TCM1::TLastCell_, Cell<DPDirectionZero, False> >::Type TResult1c;

        SEQAN_ASSERT_EQ(+TResult1a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile2, DPBandFirstPhase> TCM2;

        typedef IsSameType<TCM2::TFirstCell_, Cell<DPDirectionZero, False> >::Type TResult2a;
        typedef IsSameType<TCM2::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult2b;
        typedef IsSameType<TCM2::TLastCell_, Cell<DPDirectionLowerBand, False> >::Type TResult2c;

        SEQAN_ASSERT_EQ(+TResult2a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult2b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult2c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile2, DPBandMiddlePhase> TCM3;

        typedef IsSameType<TCM3::TFirstCell_, Cell<DPDirectionUpperBand, False> >::Type TResult3a;
        typedef IsSameType<TCM3::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult3b;
        typedef IsSameType<TCM3::TLastCell_, Cell<DPDirectionLowerBand, False> >::Type TResult3c;

        SEQAN_ASSERT_EQ(+TResult3a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile2, DPBandLastPhase> TCM4;

        typedef IsSameType<TCM4::TFirstCell_, Cell<DPDirectionUpperBand, False> >::Type TResult4a;
        typedef IsSameType<TCM4::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult4b;
        typedef IsSameType<TCM4::TLastCell_, Cell<DPDirectionAll, True> >::Type TResult4c;

        SEQAN_ASSERT_EQ(+TResult4a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult4b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult4c::VALUE, true);
    }

    {   // standard local and unbanded
        typedef SetUpColumnManager<TAligmentProfile3, DPNoBandInitPhase> TCM1;

        typedef IsSameType<TCM1::TFirstCell_, Cell<DPDirectionZero, True> >::Type TResult1a;
        typedef IsSameType<TCM1::TInnerCell_, Cell<DPDirectionZero, True> >::Type TResult1b;
        typedef IsSameType<TCM1::TLastCell_, Cell<DPDirectionZero, True> >::Type TResult1c;

        SEQAN_ASSERT_EQ(+TResult1a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile3, DPNoBandPhase> TCM3;

        typedef IsSameType<TCM3::TFirstCell_, Cell<DPDirectionZero, True> >::Type TResult3a;
        typedef IsSameType<TCM3::TInnerCell_, Cell<DPDirectionAll, True> >::Type TResult3b;
        typedef IsSameType<TCM3::TLastCell_, Cell<DPDirectionAll, True> >::Type TResult3c;

        SEQAN_ASSERT_EQ(+TResult3a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3c::VALUE, true);
    }

    {   // standard local and banded
        typedef SetUpColumnManager<TAligmentProfile3, DPBandInitPhase> TCM1;

        typedef IsSameType<TCM1::TFirstCell_, Cell<DPDirectionZero, True> >::Type TResult1a;
        typedef IsSameType<TCM1::TInnerCell_, Cell<DPDirectionZero, True> >::Type TResult1b;
        typedef IsSameType<TCM1::TLastCell_, Cell<DPDirectionZero, True> >::Type TResult1c;

        SEQAN_ASSERT_EQ(+TResult1a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile3, DPBandFirstPhase> TCM2;

        typedef IsSameType<TCM2::TFirstCell_, Cell<DPDirectionZero, True> >::Type TResult2a;
        typedef IsSameType<TCM2::TInnerCell_, Cell<DPDirectionAll, True> >::Type TResult2b;
        typedef IsSameType<TCM2::TLastCell_, Cell<DPDirectionLowerBand, True> >::Type TResult2c;

        SEQAN_ASSERT_EQ(+TResult2a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult2b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult2c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile3, DPBandMiddlePhase> TCM3;

        typedef IsSameType<TCM3::TFirstCell_, Cell<DPDirectionUpperBand, True> >::Type TResult3a;
        typedef IsSameType<TCM3::TInnerCell_, Cell<DPDirectionAll, True> >::Type TResult3b;
        typedef IsSameType<TCM3::TLastCell_, Cell<DPDirectionLowerBand, True> >::Type TResult3c;

        SEQAN_ASSERT_EQ(+TResult3a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile3, DPBandLastPhase> TCM4;

        typedef IsSameType<TCM4::TFirstCell_, Cell<DPDirectionUpperBand, True> >::Type TResult4a;
        typedef IsSameType<TCM4::TInnerCell_, Cell<DPDirectionAll, True> >::Type TResult4b;
        typedef IsSameType<TCM4::TLastCell_, Cell<DPDirectionAll, True> >::Type TResult4c;

        SEQAN_ASSERT_EQ(+TResult4a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult4b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult4c::VALUE, true);
    }

    {   // global with row max and unbanded
        typedef SetUpColumnManager<TAligmentProfile4, DPNoBandInitPhase> TCM1;

        typedef IsSameType<TCM1::TFirstCell_, Cell<DPDirectionZero, True> >::Type TResult1a;
        typedef IsSameType<TCM1::TInnerCell_, Cell<DPDirectionVertical, True> >::Type TResult1b;
        typedef IsSameType<TCM1::TLastCell_, Cell<DPDirectionVertical, True> >::Type TResult1c;

        SEQAN_ASSERT_EQ(+TResult1a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile4, DPNoBandPhase> TCM3;

        typedef IsSameType<TCM3::TFirstCell_, Cell<DPDirectionHorizontal, True> >::Type TResult3a;
        typedef IsSameType<TCM3::TInnerCell_, Cell<DPDirectionAll, True> >::Type TResult3b;
        typedef IsSameType<TCM3::TLastCell_, Cell<DPDirectionAll, True> >::Type TResult3c;

        SEQAN_ASSERT_EQ(+TResult3a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3c::VALUE, true);
    }

    {  // global with row max and banded
        typedef SetUpColumnManager<TAligmentProfile4, DPBandInitPhase> TCM1;

        typedef IsSameType<TCM1::TFirstCell_, Cell<DPDirectionZero, True> >::Type TResult1a;
        typedef IsSameType<TCM1::TInnerCell_, Cell<DPDirectionVertical, True> >::Type TResult1b;
        typedef IsSameType<TCM1::TLastCell_, Cell<DPDirectionVertical, True> >::Type TResult1c;

        SEQAN_ASSERT_EQ(+TResult1a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile4, DPBandFirstPhase> TCM2;

        typedef IsSameType<TCM2::TFirstCell_, Cell<DPDirectionHorizontal, True> >::Type TResult2a;
        typedef IsSameType<TCM2::TInnerCell_, Cell<DPDirectionAll, True> >::Type TResult2b;
        typedef IsSameType<TCM2::TLastCell_, Cell<DPDirectionLowerBand, True> >::Type TResult2c;

        SEQAN_ASSERT_EQ(+TResult2a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult2b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult2c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile4, DPBandMiddlePhase> TCM3;

        typedef IsSameType<TCM3::TFirstCell_, Cell<DPDirectionUpperBand, True> >::Type TResult3a;
        typedef IsSameType<TCM3::TInnerCell_, Cell<DPDirectionAll, True> >::Type TResult3b;
        typedef IsSameType<TCM3::TLastCell_, Cell<DPDirectionLowerBand, True> >::Type TResult3c;

        SEQAN_ASSERT_EQ(+TResult3a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile4, DPBandLastPhase> TCM4;

        typedef IsSameType<TCM4::TFirstCell_, Cell<DPDirectionUpperBand, True> >::Type TResult4a;
        typedef IsSameType<TCM4::TInnerCell_, Cell<DPDirectionAll, True> >::Type TResult4b;
        typedef IsSameType<TCM4::TLastCell_, Cell<DPDirectionAll, True> >::Type TResult4c;

        SEQAN_ASSERT_EQ(+TResult4a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult4b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult4c::VALUE, true);
    }
}

void testAlignmentManager()
{
    typedef Band<BandSwitchedOff> TBand1;
    typedef Band<BandSwitchedOn<> > TBand2;

    typedef AlignmentProfile<Global<FreeEndGaps<> >, LinearGaps, TracebackSwitchedOff > TAlignmentProfile1;

    typedef DPManager<TAlignmentProfile1, TBand1> TDPManager1;
    typedef DPManager<TAlignmentProfile1, TBand2> TDPManager2;


    TDPManager1 cm;

    SEQAN_ASSERT_EQ(cm._spanDp, 0);
    SEQAN_ASSERT_EQ(cm._spanTrace, 0);
    SEQAN_ASSERT_EQ(cm._spanSeqVBegin, 0);
    SEQAN_ASSERT_EQ(cm._spanSeqVEnd, 0);

    cm._spanDp = 3;
    cm._spanTrace = 6;
    cm._spanSeqVBegin = 12;
    cm._spanSeqVEnd = 8;

    TDPManager1 cm2(cm);
    SEQAN_ASSERT_EQ(cm2._spanDp, 3);
    SEQAN_ASSERT_EQ(cm2._spanTrace, 6);
    SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, 12);
    SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 8);

    TDPManager1 cm3(8, 4);

    SEQAN_ASSERT_EQ(cm3._spanDp, -7);
    SEQAN_ASSERT_EQ(cm3._spanTrace, 1);
    SEQAN_ASSERT_EQ(cm3._spanSeqVBegin, 0);
    SEQAN_ASSERT_EQ(cm3._spanSeqVEnd, 6);

    TDPManager2 cm4(8, 4);

    SEQAN_ASSERT_EQ(cm4._spanDp, -3);
    SEQAN_ASSERT_EQ(cm4._spanTrace, 5);
    SEQAN_ASSERT_EQ(cm4._spanSeqVBegin, -1);
    SEQAN_ASSERT_EQ(cm4._spanSeqVEnd, 2);
}

// TODO(rmaerker): remove deprecated code
void testAlignmentManagerGetColumnManager()
{
//    typedef Band<BandSwitchedOff> TBand1;
//
//    typedef AlignmentProfile<Global<FreeEndGaps<> >, LinearGaps, TracebackSwitchedOff > TAlignmentProfile1;
//    typedef DPManager<TAlignmentProfile1, TBand1> TDPManager1;
//
//    typedef IsSameType<GetColumnManager<TDPManager1, DPPhase>::Type, TDPManager1::TInitColumnManager_>::Type TResult1;
//    typedef IsSameType<GetColumnManager<TDPManager1, BandOpenColumn>::Type, TDPManager1::TBeginColumnManager_>::Type TResult2;
//    typedef IsSameType<GetColumnManager<TDPManager1, FullColumn>::Type, TDPManager1::TMiddleColumnManager_>::Type TResult3;
//    typedef IsSameType<GetColumnManager<TDPManager1, BandCloseColumn>::Type, TDPManager1::TEndColumnManager_>::Type TResult4;
//
//    SEQAN_ASSERT_EQ(+TResult1::VALUE, true);
//    SEQAN_ASSERT_EQ(+TResult2::VALUE, true);
//    SEQAN_ASSERT_EQ(+TResult3::VALUE, true);
//    SEQAN_ASSERT_EQ(+TResult4::VALUE, true);
}

void testAlignmentManagerInit()
{
    typedef Band<BandSwitchedOff> TBand1;
    typedef Band<BandSwitchedOn<> > TBand2;

    typedef AlignmentProfile<Global<FreeEndGaps<> >, LinearGaps, TracebackSwitchedOff > TAlignmentProfile1;

    typedef DPManager<TAlignmentProfile1, TBand1> TDPManager1;
    typedef DPManager<TAlignmentProfile1, TBand2> TDPManager2;

    TDPManager1 cm1;
    _init(cm1, 8, 4);

    SEQAN_ASSERT_EQ(cm1._spanDp, -7);
    SEQAN_ASSERT_EQ(cm1._spanTrace, 1);
    SEQAN_ASSERT_EQ(cm1._spanSeqVBegin, 0);
    SEQAN_ASSERT_EQ(cm1._spanSeqVEnd, 6);

    TDPManager2 cm2;
    _init(cm2, 8, 4);

    SEQAN_ASSERT_EQ(cm2._spanDp, -3);
    SEQAN_ASSERT_EQ(cm2._spanTrace, 5);
    SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, -1);
    SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 2);
}

void testAlignmentManagerUpdate()
{
    typedef Band<BandSwitchedOff> TBand1;
    typedef Band<BandSwitchedOn<> > TBand2;

    typedef AlignmentProfile<Global<>, LinearGaps, TracebackSwitchedOff > TAlignmentProfile1;

    typedef DPManager<TAlignmentProfile1, TBand1> TDPManager1;
    typedef DPManager<TAlignmentProfile1, TBand2> TDPManager2;

    {
        TDPManager1 cm1(8, 4);

        SEQAN_ASSERT_EQ(cm1._spanDp, -7);
        SEQAN_ASSERT_EQ(cm1._spanTrace, 1);
        SEQAN_ASSERT_EQ(cm1._spanSeqVBegin, 0);
        SEQAN_ASSERT_EQ(cm1._spanSeqVEnd, 6);

        update(cm1, DPBandFirstPhase());
        SEQAN_ASSERT_EQ(cm1._spanDp, -7);
        SEQAN_ASSERT_EQ(cm1._spanTrace, 1);
        SEQAN_ASSERT_EQ(cm1._spanSeqVBegin, 0);
        SEQAN_ASSERT_EQ(cm1._spanSeqVEnd, 6);


        for(unsigned int i = 1; i < 8; ++i)
        {
            update(cm1, DPNoBandPhase());
            SEQAN_ASSERT_EQ(cm1._spanDp, -7);
            SEQAN_ASSERT_EQ(cm1._spanTrace, 1);
            SEQAN_ASSERT_EQ(cm1._spanSeqVBegin, 0);
            SEQAN_ASSERT_EQ(cm1._spanSeqVEnd, 6);
        }
    }


    {
        TDPManager2 cm2(8, 4);

        SEQAN_ASSERT_EQ(cm2._spanDp, -3);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 5);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, -1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 2);

        update(cm2, DPBandInitPhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -3);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 5);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, -1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 2);

        update(cm2, DPBandFirstPhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -4);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 4);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, -1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 3);

        update(cm2, DPBandFirstPhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -5);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 3);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, -1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 4);

        update(cm2, DPBandFirstPhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -6);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 2);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, -1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 5);

        update(cm2, DPBandFirstPhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -7);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, -1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 6);

        update(cm2, DPBandMiddlePhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -7);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, 0);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 7);

        update(cm2, DPBandMiddlePhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -7);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, 1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 8);

        update(cm2, DPBandMiddlePhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -7);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, 2);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 9);

        update(cm2, DPBandMiddlePhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -7);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 1);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, 3);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 10);

        update(cm2, DPBandLastPhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -6);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 2);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, 4);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 10);

        update(cm2, DPBandLastPhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -5);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 3);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, 5);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 10);


        update(cm2, DPBandLastPhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -4);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 4);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, 6);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 10);

        update(cm2, DPBandLastPhase());
        SEQAN_ASSERT_EQ(cm2._spanDp, -3);
        SEQAN_ASSERT_EQ(cm2._spanTrace, 5);
        SEQAN_ASSERT_EQ(cm2._spanSeqVBegin, 7);
        SEQAN_ASSERT_EQ(cm2._spanSeqVEnd, 10);
    }

}

void
testAlignmentManagerBandedChain()
{
    typedef Band<BandSwitchedOff> TBand1;
    typedef Band<BandSwitchedOn<> > TBand2;

    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOff > TAligmentProfile;

    { // standard global and unbanded
        typedef SetUpColumnManager<TAligmentProfile, DPNoBandInitPhase> TCM1;

        typedef IsSameType<TCM1::TFirstCell_, Cell<DPDirectionAdopt, False> >::Type TResult1a;
        typedef IsSameType<TCM1::TInnerCell_, Cell<DPDirectionAdopt, False> >::Type TResult1b;
        typedef IsSameType<TCM1::TLastCell_, Cell<DPDirectionAdopt, True> >::Type TResult1c;

        SEQAN_ASSERT_EQ(+TResult1a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile, DPNoBandPhase> TCM3;

        typedef IsSameType<TCM3::TFirstCell_, Cell<DPDirectionAdopt, False> >::Type TResult3a;
        typedef IsSameType<TCM3::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult3b;
        typedef IsSameType<TCM3::TLastCell_, Cell<DPDirectionAll, True> >::Type TResult3c;

        SEQAN_ASSERT_EQ(+TResult3a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3c::VALUE, true);
    }

    { // standard global and unbanded
        typedef SetUpColumnManager<TAligmentProfile, DPBandInitPhase> TCM1;

        typedef IsSameType<TCM1::TFirstCell_, Cell<DPDirectionAdopt, False> >::Type TResult1a;
        typedef IsSameType<TCM1::TInnerCell_, Cell<DPDirectionAdopt, False> >::Type TResult1b;
        typedef IsSameType<TCM1::TLastCell_, Cell<DPDirectionAdopt, False> >::Type TResult1c;

        SEQAN_ASSERT_EQ(+TResult1a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult1c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile, DPBandFirstPhase> TCM2;

        typedef IsSameType<TCM2::TFirstCell_, Cell<DPDirectionAdopt, False> >::Type TResult2a;
        typedef IsSameType<TCM2::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult2b;
        typedef IsSameType<TCM2::TLastCell_, Cell<DPDirectionLowerBand, False> >::Type TResult2c;

        SEQAN_ASSERT_EQ(+TResult2a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult2b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult2c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile, DPBandMiddlePhase> TCM3;

        typedef IsSameType<TCM3::TFirstCell_, Cell<DPDirectionUpperBand, False> >::Type TResult3a;
        typedef IsSameType<TCM3::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult3b;
        typedef IsSameType<TCM3::TLastCell_, Cell<DPDirectionLowerBand, False> >::Type TResult3c;

        SEQAN_ASSERT_EQ(+TResult3a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult3c::VALUE, true);

        typedef SetUpColumnManager<TAligmentProfile, DPBandLastPhase> TCM4;

        typedef IsSameType<TCM4::TFirstCell_, Cell<DPDirectionUpperBand, False> >::Type TResult4a;
        typedef IsSameType<TCM4::TInnerCell_, Cell<DPDirectionAll, False> >::Type TResult4b;
        typedef IsSameType<TCM4::TLastCell_, Cell<DPDirectionAll, True> >::Type TResult4c;

        SEQAN_ASSERT_EQ(+TResult4a::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult4b::VALUE, true);
        SEQAN_ASSERT_EQ(+TResult4c::VALUE, true);
    }
}

SEQAN_DEFINE_TEST(test_alignment_dp_manager_get_dp_direction)
{
    testAlignmentManagerGetDpDirection();
}

SEQAN_DEFINE_TEST(test_alignment_dp_manager_is_to_track)
{
    testAlignmentManagerIsToTrack();
}

SEQAN_DEFINE_TEST(test_alignment_dp_manager_setup_column_manager)
{
    testAlignmentManagerSetUpColumnManager();
}

SEQAN_DEFINE_TEST(test_alignment_dp_manager)
{
    testAlignmentManager();
}

SEQAN_DEFINE_TEST(test_alignment_dp_manager_get_column_manager)
{
    testAlignmentManagerGetColumnManager();
}


SEQAN_DEFINE_TEST(test_alignment_dp_manager_init)
{
    testAlignmentManagerInit();
}

SEQAN_DEFINE_TEST(test_alignment_dp_manager_update)
{
    testAlignmentManagerUpdate();
}

SEQAN_DEFINE_TEST(test_alignment_dp_manager_banded_chain)
{
    testAlignmentManagerBandedChain();
}

#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_MANAGER_H_
