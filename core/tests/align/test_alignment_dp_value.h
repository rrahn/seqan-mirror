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

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_VALUE_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_VALUE_H_

using namespace seqan;

#include <seqan/basic.h>

#include <seqan/align/alignment_base.h>
#include <seqan/align/alignment_dp_value.h>


void testAlignmentDpValueLinearConstructor()
{
    typedef DPValue<int, LinearGaps, DPValueConfigDefault> TDPValue;
    { // Default Ctor
        TDPValue dpValue;

        SEQAN_ASSERT_EQ(dpValue._score, 0);
    }

    { // Copy Ctor
        TDPValue dpValue;
        dpValue._score = -10;
        TDPValue dpValue2(dpValue);
        SEQAN_ASSERT_EQ(dpValue2._score, -10);
    }

    { // Ctor with score
        TDPValue dpValue(8);
        SEQAN_ASSERT_EQ(dpValue._score, 8);
    }
}

void testAlignmentDpValueLinearForbiddenConstructor()
{
    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue;
    { // Default Ctor
        TDPValue dpValue;

        SEQAN_ASSERT_EQ(dpValue._score, 0);
        SEQAN_ASSERT_EQ(dpValue._forbidden, false);
    }

    { // Copy Ctor
        TDPValue dpValue;
        dpValue._score = -10;
        TDPValue dpValue2(dpValue);
        SEQAN_ASSERT_EQ(dpValue2._score, -10);
        SEQAN_ASSERT_EQ(dpValue2._forbidden, false);

        dpValue._forbidden = true;
        TDPValue dpValue3 = dpValue;
        SEQAN_ASSERT_EQ(dpValue3._score, -10);
        SEQAN_ASSERT_EQ(dpValue3._forbidden, true);
    }

    { // Ctor with score
        TDPValue dpValue(8);
        SEQAN_ASSERT_EQ(dpValue._score, 8);
        SEQAN_ASSERT_EQ(dpValue._forbidden, false);
    }

    { // Ctor with forbidden flag
        TDPValue dpValue(true);
        SEQAN_ASSERT_EQ(dpValue._score, 0);
        SEQAN_ASSERT_EQ(dpValue._forbidden, true);

        TDPValue dpValue2(false);
        SEQAN_ASSERT_EQ(dpValue2._score, 0);
        SEQAN_ASSERT_EQ(dpValue2._forbidden, false);
    }
}

void testAlignmentDpValueAffineConstructor()
{
    typedef DPValue<int, AffineGaps, DPValueConfigDefault> TDPValue;
    { // Default Ctor
        TDPValue dpValue;

        SEQAN_ASSERT_EQ(dpValue._score, 0);
        SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
        SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
        SEQAN_ASSERT_EQ(dpValue._dirForbidden, (+ForbiddenDirection::NONE));
    }

    { // Copy Ctor
        TDPValue dpValue;
        dpValue._score = -10;
        dpValue._scoreHorizontal = 6;
        dpValue._scoreVertical = 3;
        dpValue._dirForbidden = +ForbiddenDirection::HORIZONTAL;

        TDPValue dpValue2(dpValue);
        SEQAN_ASSERT_EQ(dpValue2._score, -10);
        SEQAN_ASSERT_EQ(dpValue2._scoreHorizontal, 6);
        SEQAN_ASSERT_EQ(dpValue2._scoreVertical, 3);
        SEQAN_ASSERT_EQ(dpValue2._dirForbidden, (+ForbiddenDirection::HORIZONTAL));
    }

    { // Ctor with score
        TDPValue dpValue(8);

        SEQAN_ASSERT_EQ(dpValue._score, 8);
        SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
        SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
        SEQAN_ASSERT_EQ(dpValue._dirForbidden, (+ForbiddenDirection::NONE));
    }
}

void testAlignmentDpValueAffineForbiddenConstructor()
{
    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue;
    { // Default Ctor
        TDPValue dpValue;

        SEQAN_ASSERT_EQ(dpValue._score, 0);
        SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
        SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
        SEQAN_ASSERT_EQ(dpValue._dirForbidden, (+ForbiddenDirection::NONE));
        SEQAN_ASSERT_EQ(dpValue._forbidden, false);
    }

    { // Copy Ctor
        TDPValue dpValue;
        dpValue._score = -10;
        dpValue._scoreHorizontal = 6;
        dpValue._scoreVertical = 3;
        dpValue._dirForbidden = +ForbiddenDirection::HORIZONTAL;

        TDPValue dpValue2(dpValue);
        SEQAN_ASSERT_EQ(dpValue2._score, -10);
        SEQAN_ASSERT_EQ(dpValue2._scoreHorizontal, 6);
        SEQAN_ASSERT_EQ(dpValue2._scoreVertical, 3);
        SEQAN_ASSERT_EQ(dpValue2._dirForbidden, (+ForbiddenDirection::HORIZONTAL));
        SEQAN_ASSERT_EQ(dpValue2._forbidden, false);

        dpValue._forbidden = true;
        TDPValue dpValue3 = dpValue;
        SEQAN_ASSERT_EQ(dpValue3._score, -10);
        SEQAN_ASSERT_EQ(dpValue3._scoreHorizontal, 6);
        SEQAN_ASSERT_EQ(dpValue3._scoreVertical, 3);
        SEQAN_ASSERT_EQ(dpValue3._dirForbidden, (+ForbiddenDirection::HORIZONTAL));
        SEQAN_ASSERT_EQ(dpValue3._forbidden, true);
    }

    { // Ctor with score
        TDPValue dpValue(8);

        SEQAN_ASSERT_EQ(dpValue._score, 8);
        SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
        SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
        SEQAN_ASSERT_EQ(dpValue._dirForbidden, (+ForbiddenDirection::NONE));
        SEQAN_ASSERT_EQ(dpValue._forbidden, false);
    }

    { // Ctor with forbidden flag
        TDPValue dpValue(true);
        SEQAN_ASSERT_EQ(dpValue._score, 0);
        SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
        SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
        SEQAN_ASSERT_EQ(dpValue._dirForbidden, (+ForbiddenDirection::NONE));
        SEQAN_ASSERT_EQ(dpValue._forbidden, true);

        TDPValue dpValue2(false);
        SEQAN_ASSERT_EQ(dpValue2._score, 0);
        SEQAN_ASSERT_EQ(dpValue2._scoreHorizontal, 0);
        SEQAN_ASSERT_EQ(dpValue2._scoreVertical, 0);
        SEQAN_ASSERT_EQ(dpValue2._dirForbidden, (+ForbiddenDirection::NONE));
        SEQAN_ASSERT_EQ(dpValue2._forbidden, false);
    }
}

template <typename TDPValue>
void testAlignmentDpValueGetScore(TDPValue const & dpValue)
{
    SEQAN_ASSERT_EQ(getScore(dpValue), dpValue._score);
}


template <typename TDPValue, typename TScoreValue>
void testAlignmentDpValueSetScore(TDPValue & dpValue, TScoreValue const & score)
{
    SEQAN_ASSERT_NEQ(dpValue._score, score);
    setScore(dpValue, score);
    SEQAN_ASSERT_EQ(dpValue._score, score);
}

template <typename TScoreValue, typename TDPValueSpec>
void testAlignmentDpValueIsForbidden(DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> const & dpValue)
{
    SEQAN_ASSERT_EQ(isForbidden(dpValue), dpValue._forbidden);
}

template <typename TScoreValue, typename TDPValueSpec>
void testAlignmentDpValueIsForbidden(DPValue<TScoreValue, TDPValueSpec, DPValueConfigDefault> const & dpValue)
{
    SEQAN_ASSERT_EQ(isForbidden(dpValue), false);
}

template <typename TScoreValue, typename TDPValueSpec>
void testAlignmentDpValueSetForbidden(DPValue<TScoreValue, TDPValueSpec, DPValueConfigDefault> const &, bool const & forbidden)
{
    DPValue<TScoreValue, TDPValueSpec, DPValueConfigDefault> dpValue;
    setForbidden(dpValue, forbidden);
    SEQAN_ASSERT_EQ(isForbidden(dpValue), false);
}

template <typename TScoreValue, typename TDPValueSpec>
void testAlignmentDpValueSetForbidden(DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> const &, bool const & forbidden)
{
    DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> dpValue;
    SEQAN_ASSERT_EQ(dpValue._forbidden, false);
    setForbidden(dpValue, forbidden);
    SEQAN_ASSERT_EQ(dpValue._forbidden, forbidden);
}

SEQAN_DEFINE_TEST(test_alignment_dp_value_constructor)
{
    testAlignmentDpValueLinearConstructor();
    testAlignmentDpValueLinearForbiddenConstructor();
    testAlignmentDpValueAffineConstructor();
    testAlignmentDpValueAffineForbiddenConstructor();
}

SEQAN_DEFINE_TEST(test_alignment_dp_value_get_score)
{
    typedef DPValue<int, LinearGaps> TDPValue1;
    testAlignmentDpValueGetScore(TDPValue1());
    testAlignmentDpValueGetScore(TDPValue1(2));
    testAlignmentDpValueGetScore(TDPValue1(-4));
    typedef DPValue<int, LinearGaps> TDPValue2;
    testAlignmentDpValueGetScore(TDPValue2());
    testAlignmentDpValueGetScore(TDPValue2(2));
    testAlignmentDpValueGetScore(TDPValue2(-4));

    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue3;
    testAlignmentDpValueGetScore(TDPValue3());
    TDPValue3 dpValue31(2);
    dpValue31._forbidden = true;
    testAlignmentDpValueGetScore(dpValue31);
    TDPValue3 dpValue32(-4);
    dpValue32._forbidden = false;
    testAlignmentDpValueGetScore(dpValue32);

    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue4;
    testAlignmentDpValueGetScore(TDPValue4());
    TDPValue4 dpValue41(2);
    dpValue41._scoreHorizontal = 12;
    dpValue41._scoreVertical = 12;
    dpValue41._forbidden = true;
    testAlignmentDpValueGetScore(dpValue41);
}

SEQAN_DEFINE_TEST(test_alignment_dp_value_set_score)
{
    typedef DPValue<int, LinearGaps> TDPValue1;
    TDPValue1 dpValue11(3);
    testAlignmentDpValueSetScore(dpValue11, 2);
    TDPValue1 dpValue12(3);
    testAlignmentDpValueSetScore(dpValue12, -4);

    typedef DPValue<int, AffineGaps> TDPValue2;
    TDPValue2 dpValue21(3);
    testAlignmentDpValueSetScore(dpValue21, 2);
    TDPValue2 dpValue22(3);
    testAlignmentDpValueSetScore(dpValue22, -4);

    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue3;
    TDPValue3 dpValue31(3);
    dpValue31._forbidden = false;
    TDPValue3 dpValue32(3);
    dpValue32._forbidden = true;
    testAlignmentDpValueSetScore(dpValue32, -4);

    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue4;
    TDPValue4 dpValue41(3);
    dpValue41._forbidden = false;
    testAlignmentDpValueSetScore(dpValue41, 2);
    TDPValue4 dpValue42(3);
    dpValue42._forbidden = true;
    testAlignmentDpValueSetScore(dpValue42, -4);

}

SEQAN_DEFINE_TEST(test_alignment_dp_value_get_origin)
{
    typedef DPValue<int, LinearGaps> TDPValue1;

    testAlignmentDpValueIsForbidden(TDPValue1());
    TDPValue1 dpValue11(3);
    testAlignmentDpValueIsForbidden(dpValue11);
    TDPValue1 dpValue12(3);
    testAlignmentDpValueIsForbidden(dpValue12);

    typedef DPValue<int, AffineGaps> TDPValue2;
    testAlignmentDpValueIsForbidden(TDPValue2());
    TDPValue2 dpValue21(3);
    testAlignmentDpValueIsForbidden(dpValue21);
    TDPValue2 dpValue22(3);
    testAlignmentDpValueIsForbidden(dpValue22);


    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue3;

    testAlignmentDpValueIsForbidden(TDPValue3());
    TDPValue3 dpValue31(3);
    dpValue31._forbidden = false;
    testAlignmentDpValueIsForbidden(dpValue31);
    TDPValue3 dpValue32(3);
    dpValue32._forbidden = true;
    testAlignmentDpValueIsForbidden(dpValue32);

    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue4;
    testAlignmentDpValueIsForbidden(TDPValue4());
    TDPValue4 dpValue41(3);
    dpValue41._forbidden = false;
    testAlignmentDpValueIsForbidden(dpValue41);
    TDPValue4 dpValue42(3);
    dpValue42._forbidden = true;
    testAlignmentDpValueIsForbidden(dpValue42);
}

SEQAN_DEFINE_TEST(test_alignment_dp_value_set_origin)
{
    testAlignmentDpValueSetForbidden(DPValue<int, LinearGaps>(), true);
    testAlignmentDpValueSetForbidden(DPValue<int, LinearGaps>(), false);
    testAlignmentDpValueSetForbidden(DPValue<int, AffineGaps>(), true);
    testAlignmentDpValueSetForbidden(DPValue<int, AffineGaps>(), false);
    testAlignmentDpValueSetForbidden(DPValue<int, LinearGaps, DPValueConfigForbidden>(), true);
    testAlignmentDpValueSetForbidden(DPValue<int, LinearGaps, DPValueConfigForbidden>(), false);
    testAlignmentDpValueSetForbidden(DPValue<int, AffineGaps, DPValueConfigForbidden>(), true);
    testAlignmentDpValueSetForbidden(DPValue<int, AffineGaps, DPValueConfigForbidden>(), false);
}

#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_VALUE_H_
