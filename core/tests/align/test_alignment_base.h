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

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_BASE_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_BASE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/align/alignment_base.h>

using namespace seqan;

template <typename TCoreSpec, typename TGapSpec, typename TTraceSpec>
inline void testalignmentBaseGetAlignmentSpec(AlignmentProfile<TCoreSpec, TGapSpec, TTraceSpec> const &)
{
    typedef AlignmentProfile<TCoreSpec, TGapSpec, TTraceSpec> TAlignmentProfile;

    typedef typename IsSameType<TCoreSpec, typename GetAlignmentSpec<TAlignmentProfile>::Type >::Type TResult;
    SEQAN_ASSERT_EQ(+TResult::VALUE, true);
}

template <typename TCoreSpec, typename TGapSpec, typename TTraceSpec>
inline void testalignmentBaseGetGapSpec(AlignmentProfile<TCoreSpec, TGapSpec, TTraceSpec> const &)
{
    typedef AlignmentProfile<TCoreSpec, TGapSpec, TTraceSpec>  TAlignmentProfile;

    typedef typename IsSameType<TGapSpec, typename GetGapSpec<TAlignmentProfile>::Type >::Type TResult1;

    SEQAN_ASSERT_EQ(+TResult1::VALUE, true);
}

template <typename TCoreSpec, typename TGapSpec, typename TTraceSpec>
inline void testalignmentBaseIsTracebackOn(AlignmentProfile<TCoreSpec, TGapSpec, TTraceSpec> const &)
{
    typedef AlignmentProfile<TCoreSpec, TGapSpec, TTraceSpec>  TAlignmentProfile;

    typedef typename IsSameType<TTraceSpec, TracebackSwitchedOn >::Type TCompare;
    if (TCompare::VALUE)
    {
        typedef typename IsTracebackOn<TAlignmentProfile>::Type TIsOn;
        typedef typename IsSameType<TIsOn, True>::Type TResult;
        SEQAN_ASSERT_EQ(+TResult::VALUE, true);
    }
    else
    {
        typedef typename IsTracebackOn<TAlignmentProfile>::Type TIsOn;
        typedef typename IsSameType<TIsOn, True>::Type TResult;
        SEQAN_ASSERT_EQ(+TResult::VALUE, true);
    }
}

template <typename TCoreSpec, typename TGapSpec, typename TTraceSpec>
inline void testalignmentBaseIsLocal(AlignmentProfile<TCoreSpec, TGapSpec, TTraceSpec> const &)
{
    typedef AlignmentProfile<TCoreSpec, TGapSpec, TTraceSpec > TAlignmentProfile;
    typedef typename Spec<TCoreSpec>::Type TSpec_;

    typedef typename IsSameType<typename IsLocal<TAlignmentProfile>::Type,
            typename IsSameType<Local<TSpec_>, TCoreSpec>::Type >::Type TResult1;

    SEQAN_ASSERT_EQ(+TResult1::VALUE, true);
}

template <typename TCoreSpec, typename TGapSpec, typename TTraceSpec>
inline void testalignmentBaseIsGlobal(AlignmentProfile<TCoreSpec, TGapSpec, TTraceSpec> const &)
{
    typedef AlignmentProfile<TCoreSpec, TGapSpec, TTraceSpec> TAlignmentProfile;
    typedef typename Spec<TCoreSpec>::Type TSpec_;

    typedef typename IsSameType<typename IsGlobal<TAlignmentProfile>::Type,
        typename IsSameType<Global<TSpec_>, TCoreSpec>::Type >::Type TResult1;

    SEQAN_ASSERT_EQ(+TResult1::VALUE, true);
}

template<typename TAlignmentCore>
void testalignmentBaseIsFreeEndGaps()
{
    {
        typedef AlignmentProfile<Global<FreeEndGaps<> >, LinearGaps, TracebackSwitchedOff > TAlignmentProfile;

        bool test1 = IsFreeEndGap<TAlignmentProfile, FirstRow>::VALUE;
        bool test2 = IsFreeEndGap<TAlignmentProfile, FirstColumn>::VALUE;
        bool test3 = IsFreeEndGap<TAlignmentProfile, LastRow>::VALUE;
        bool test4 = IsFreeEndGap<TAlignmentProfile, LastColumn>::VALUE;

        SEQAN_ASSERT_NOT(test1);
        SEQAN_ASSERT_NOT(test2);
        SEQAN_ASSERT_NOT(test3);
        SEQAN_ASSERT_NOT(test4);
    }

    {
        typedef AlignmentProfile<Global<FreeEndGaps<True> >, LinearGaps, TracebackSwitchedOff > TAlignmentProfile;

        bool test1 = IsFreeEndGap<TAlignmentProfile, FirstRow>::VALUE;
        bool test2 = IsFreeEndGap<TAlignmentProfile, FirstColumn>::VALUE;
        bool test3 = IsFreeEndGap<TAlignmentProfile, LastRow>::VALUE;
        bool test4 = IsFreeEndGap<TAlignmentProfile, LastColumn>::VALUE;

        SEQAN_ASSERT(test1);
        SEQAN_ASSERT_NOT(test2);
        SEQAN_ASSERT_NOT(test3);
        SEQAN_ASSERT_NOT(test4);
    }

    {
        typedef AlignmentProfile<Global<FreeEndGaps<True, True> >, LinearGaps, TracebackSwitchedOff > TAlignmentProfile;

        bool test1 = IsFreeEndGap<TAlignmentProfile, FirstRow>::VALUE;
        bool test2 = IsFreeEndGap<TAlignmentProfile, FirstColumn>::VALUE;
        bool test3 = IsFreeEndGap<TAlignmentProfile, LastRow>::VALUE;
        bool test4 = IsFreeEndGap<TAlignmentProfile, LastColumn>::VALUE;

        SEQAN_ASSERT(test1);
        SEQAN_ASSERT(test2);
        SEQAN_ASSERT_NOT(test3);
        SEQAN_ASSERT_NOT(test4);
    }

    {
        typedef AlignmentProfile<Global<FreeEndGaps<True, True, True> >, LinearGaps, TracebackSwitchedOff > TAlignmentProfile;

        bool test1 = IsFreeEndGap<TAlignmentProfile, FirstRow>::VALUE;
        bool test2 = IsFreeEndGap<TAlignmentProfile, FirstColumn>::VALUE;
        bool test3 = IsFreeEndGap<TAlignmentProfile, LastRow>::VALUE;
        bool test4 = IsFreeEndGap<TAlignmentProfile, LastColumn>::VALUE;

        SEQAN_ASSERT(test1);
        SEQAN_ASSERT(test2);
        SEQAN_ASSERT(test3);
        SEQAN_ASSERT_NOT(test4);
    }

    {
        typedef AlignmentProfile<Global<FreeEndGaps<True, True, True, True> >, LinearGaps, TracebackSwitchedOff > TAlignmentProfile;

        bool test1 = IsFreeEndGap<TAlignmentProfile, FirstRow>::VALUE;
        bool test2 = IsFreeEndGap<TAlignmentProfile, FirstColumn>::VALUE;
        bool test3 = IsFreeEndGap<TAlignmentProfile, LastRow>::VALUE;
        bool test4 = IsFreeEndGap<TAlignmentProfile, LastColumn>::VALUE;

        SEQAN_ASSERT(test1);
        SEQAN_ASSERT(test2);
        SEQAN_ASSERT(test3);
        SEQAN_ASSERT(test4);
    }

    {
        typedef AlignmentProfile<Global<FreeEndGaps<False, True, False, True> >, LinearGaps, TracebackSwitchedOff > TAlignmentProfile;

        bool test1 = IsFreeEndGap<TAlignmentProfile, FirstRow>::VALUE;
        bool test2 = IsFreeEndGap<TAlignmentProfile, FirstColumn>::VALUE;
        bool test3 = IsFreeEndGap<TAlignmentProfile, LastRow>::VALUE;
        bool test4 = IsFreeEndGap<TAlignmentProfile, LastColumn>::VALUE;

        SEQAN_ASSERT_NOT(test1);
        SEQAN_ASSERT(test2);
        SEQAN_ASSERT_NOT(test3);
        SEQAN_ASSERT(test4);
    }
}

SEQAN_DEFINE_TEST(test_alignment_base_get_alignment_spec)
{
    typedef AlignmentProfile<Global<FreeEndGaps<True, True, False, True> >, LinearGaps, TracebackSwitchedOn> TAlignmentProfile1;
    typedef AlignmentProfile<Local<WatermanEggert>, LinearGaps, TracebackSwitchedOn > TAlignmentProfile2;
    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile3;
    typedef AlignmentProfile<Local<SplitBreakpoint>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile4;

    testalignmentBaseGetAlignmentSpec(TAlignmentProfile1());
    testalignmentBaseGetAlignmentSpec(TAlignmentProfile2());
    testalignmentBaseGetAlignmentSpec(TAlignmentProfile3());
    testalignmentBaseGetAlignmentSpec(TAlignmentProfile4());
}

SEQAN_DEFINE_TEST(test_alignment_base_get_gap_spec)
{
    typedef AlignmentProfile<Global<FreeEndGaps<True, True, False, True> >, LinearGaps, TracebackSwitchedOn > TAlignmentProfile1;
    typedef AlignmentProfile<Local<WatermanEggert>, LinearGaps, TracebackSwitchedOn > TAlignmentProfile2;
    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile3;
    typedef AlignmentProfile<Local<SplitBreakpoint>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile4;

    testalignmentBaseGetGapSpec(TAlignmentProfile1());
    testalignmentBaseGetGapSpec(TAlignmentProfile2());
    testalignmentBaseGetGapSpec(TAlignmentProfile3());
    testalignmentBaseGetGapSpec(TAlignmentProfile4());
}

SEQAN_DEFINE_TEST(test_alignment_base_is_traceback_on)
{
    typedef AlignmentProfile<Global<FreeEndGaps<True, True, False, True> >, LinearGaps, TracebackSwitchedOn > TAlignmentProfile1;
    typedef AlignmentProfile<Local<WatermanEggert>, LinearGaps, TracebackSwitchedOn > TAlignmentProfile2;
    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile3;
    typedef AlignmentProfile<Local<SplitBreakpoint>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile4;

    testalignmentBaseIsTracebackOn(TAlignmentProfile1());
    testalignmentBaseIsTracebackOn(TAlignmentProfile2());
    testalignmentBaseIsTracebackOn(TAlignmentProfile3());
    testalignmentBaseIsTracebackOn(TAlignmentProfile4());
}

SEQAN_DEFINE_TEST(test_alignment_base_is_free_end_gap)
{
    testalignmentBaseIsFreeEndGaps<AlignmentProfile<Global<> > >();
    testalignmentBaseIsFreeEndGaps<AlignmentProfile<Local<>, AffineGaps> >();
    testalignmentBaseIsFreeEndGaps<AlignmentProfile<Global<>, AffineGaps, TracebackSwitchedOn > >();
    testalignmentBaseIsFreeEndGaps<AlignmentProfile<Local<>, AffineGaps, TracebackSwitchedOff > >();
}

SEQAN_DEFINE_TEST(test_alignment_base_is_local)
{
    typedef AlignmentProfile<Global<FreeEndGaps<True, True, False, True> >, LinearGaps, TracebackSwitchedOn > TAlignmentProfile1;
    typedef AlignmentProfile<Local<WatermanEggert>, LinearGaps, TracebackSwitchedOn > TAlignmentProfile2;
    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile3;
    typedef AlignmentProfile<Local<SplitBreakpoint>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile4;

    testalignmentBaseIsLocal(TAlignmentProfile1());
    testalignmentBaseIsLocal(TAlignmentProfile2());
    testalignmentBaseIsLocal(TAlignmentProfile3());
    testalignmentBaseIsLocal(TAlignmentProfile4());
}

SEQAN_DEFINE_TEST(test_alignment_base_is_global)
{
    typedef AlignmentProfile<Global<FreeEndGaps<True, True, False, True> >, LinearGaps, TracebackSwitchedOn > TAlignmentProfile1;
    typedef AlignmentProfile<Local<WatermanEggert>, LinearGaps, TracebackSwitchedOn > TAlignmentProfile2;
    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile3;
    typedef AlignmentProfile<Local<SplitBreakpoint>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile4;

    testalignmentBaseIsGlobal(TAlignmentProfile1());
    testalignmentBaseIsGlobal(TAlignmentProfile2());
    testalignmentBaseIsGlobal(TAlignmentProfile3());
    testalignmentBaseIsGlobal(TAlignmentProfile4());
}

#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_BASE_H_
