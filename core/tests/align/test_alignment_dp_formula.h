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

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_FORMULA_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_FORMULA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/score.h>

#include <seqan/align/alignment_base.h>
#include <seqan/align/alignment_dp_band.h>
#include <seqan/align/alignment_dp_value.h>
#include <seqan/align/alignment_dp_formula.h>


void testAlignmentDpFormulaComputeScoreDiagonalLinear()
{
    using namespace seqan;

    Score<int, Simple> score(2,-2,-4);
    typedef DPValue<int, LinearGaps, DPValueConfigDefault> TDPValue;

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionDiagonal()), +TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 0);
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionDiagonal()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, -4);
}


void testAlignmentDpFormulaComputeScoreDiagonalLinearForbidden()
{
    using namespace seqan;

    Score<int, Simple> score(2,-2,-4);

    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue;
    TDPValue dpValue;
    TDPValue dpValue1(true);
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionDiagonal()), +TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
}

void testAlignmentDpFormulaComputeScoreDiagonalAffine()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigDefault> TDPValue;

    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionDiagonal()), +TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionDiagonal()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, -4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
}

void testAlignmentDpFormulaComputeScoreDiagonalAffineForbidden()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue;
    Score<int, Simple> score(2,-2,-4);


    TDPValue dpValue;
    TDPValue dpValue1(true);
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    dpValue1._scoreHorizontal = 2;
    dpValue1._scoreVertical = 2;

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionDiagonal()), +TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
}

void testAlignmentDpFormulaComputeScoreVerticalLinear()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigDefault> TDPValue;

    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionVertical()),+TraceBitMask::VERTICAL);
    SEQAN_ASSERT_EQ(dpValue._score, 2);
}

void testAlignmentDpFormulaComputeScoreVerticalLinearForbidden()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue;

    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3(true);

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionVertical()),+TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
}

void testAlignmentDpFormulaComputeScoreVerticalAffine()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigDefault> TDPValue;

    Score<int, Simple> score(2,-2,-4,-6);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);    //  6, 0, 0
    setScore(dpValue1, -2);  // -2, 0, 0
    setScore(dpValue2, 6);   //  6, 0, 0
    setScore(dpValue3, 0);   //  4, 0, 0

    dpValue3._dirForbidden = (+ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL);

    SEQAN_ASSERT_EQ(dpValue3._score, 0);
    SEQAN_ASSERT_EQ(dpValue3._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue3._scoreVertical, 0);
    SEQAN_ASSERT_EQ(dpValue3._dirForbidden, +ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionVertical()),+TraceBitMask::VERTICAL_OPEN);
    SEQAN_ASSERT_EQ(dpValue._score, -6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, -6);
    SEQAN_ASSERT_EQ(dpValue._dirForbidden, +ForbiddenDirection::HORIZONTAL);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue, 'A', 'C', score, DPDirectionVertical()),+TraceBitMask::VERTICAL);
    SEQAN_ASSERT_EQ(dpValue._score, -10);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, -10);
    SEQAN_ASSERT_EQ(dpValue._dirForbidden, +ForbiddenDirection::HORIZONTAL);
}

void testAlignmentDpFormulaComputeScoreVerticalAffineForbidden()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue;

    Score<int, Simple> score(2,-2,-4,-6);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3(true);

    setScore(dpValue, 6);    //  6, 0, 0
    setScore(dpValue1, -2);  // -2, 0, 0
    setScore(dpValue2, 6);   //  6, 0, 0
    setScore(dpValue3, 0);   //  4, 0, 0

    dpValue3._dirForbidden = (+ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL);

    SEQAN_ASSERT_EQ(dpValue3._score, 0);
    SEQAN_ASSERT_EQ(dpValue3._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue3._scoreVertical, 0);
    SEQAN_ASSERT_EQ(dpValue3._dirForbidden, +ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionVertical()),+TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(dpValue._dirForbidden, +ForbiddenDirection::HORIZONTAL);
}

void testAlignmentDpFormulaComputeScoreHorizontalLinear()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigDefault> TDPValue;
    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionHorizontal()),+TraceBitMask::HORIZONTAL);
    SEQAN_ASSERT_EQ(dpValue._score, 2);
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue, dpValue3, 'A', 'C', score, DPDirectionHorizontal()),+TraceBitMask::HORIZONTAL);
    SEQAN_ASSERT_EQ(dpValue._score, -2);
}

void testAlignmentDpFormulaComputeScoreHorizontalLinearForbidden()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue;
    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2(true);
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionHorizontal()),+TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
}

void testAlignmentDpFormulaComputeScoreHorizontalAffine()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigDefault> TDPValue;
    Score<int, Simple> score(2,-2,-4,-6);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue3, 6);


    dpValue2._dirForbidden = (+ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL);

    SEQAN_ASSERT_EQ(dpValue2._score, 0);
    SEQAN_ASSERT_EQ(dpValue2._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue2._scoreVertical, 0);
    SEQAN_ASSERT_EQ(dpValue2._dirForbidden, +ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionHorizontal()),+TraceBitMask::HORIZONTAL);
    SEQAN_ASSERT_EQ(dpValue._score, -6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, -6);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
    SEQAN_ASSERT_EQ(dpValue._dirForbidden, +ForbiddenDirection::VERTICAL);

    dpValue2 = dpValue;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionHorizontal()),+TraceBitMask::HORIZONTAL);
    SEQAN_ASSERT_EQ(dpValue._score, -10);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, -10);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
    SEQAN_ASSERT_EQ(dpValue._dirForbidden, +ForbiddenDirection::VERTICAL);
}

void testAlignmentDpFormulaComputeScoreHorizontalAffineForbidden()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue;

    Score<int, Simple> score(2,-2,-4,-6);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2(true);
    TDPValue dpValue3;

    setScore(dpValue, 6);    //  6, 0, 0
    setScore(dpValue1, -2);  // -2, 0, 0
    setScore(dpValue3, 0);   //  4, 0, 0

    dpValue2._dirForbidden = (+ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL);

    SEQAN_ASSERT_EQ(dpValue2._score, 0);
    SEQAN_ASSERT_EQ(dpValue2._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue2._scoreVertical, 0);
    SEQAN_ASSERT_EQ(dpValue2._dirForbidden, +ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionHorizontal()),+TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
    SEQAN_ASSERT_EQ(dpValue._dirForbidden, +ForbiddenDirection::VERTICAL);
}

void testAlignmentDpFormulaComputeScoreNoneLinear()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigDefault> TDPValue;

    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionZero()),+TraceBitMask::NONE);
    SEQAN_ASSERT_EQ(dpValue._score, 0);
}

void testAlignmentDpFormulaComputeScoreNoneLinearForbidden()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue;

    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1(true);
    TDPValue dpValue2(true);
    TDPValue dpValue3(true);

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionZero()),+TraceBitMask::NONE);
    SEQAN_ASSERT_EQ(dpValue._score, 0);
}

void testAlignmentDpFormulaComputeScoreNoneAffine()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigDefault> TDPValue;

    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);

    dpValue._scoreHorizontal = 10;
    dpValue._scoreVertical = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionZero()),+TraceBitMask::NONE);
    SEQAN_ASSERT_EQ(dpValue._score, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
    SEQAN_ASSERT_EQ(dpValue._dirForbidden, (+ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL));
}

void testAlignmentDpFormulaComputeScoreNoneAffineForbidden()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue;

    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, -2);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);

    dpValue._scoreHorizontal = 10;
    dpValue._scoreVertical = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionZero()),+TraceBitMask::NONE);
    SEQAN_ASSERT_EQ(dpValue._score, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
    SEQAN_ASSERT_EQ(dpValue._dirForbidden, (+ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL));

    dpValue1._forbidden = true;
    dpValue2._forbidden = true;
    dpValue3._forbidden = true;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionZero()),+TraceBitMask::NONE);
    SEQAN_ASSERT_EQ(dpValue._score, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
    SEQAN_ASSERT_EQ(dpValue._dirForbidden, (+ForbiddenDirection::HORIZONTAL | +ForbiddenDirection::VERTICAL));
}

void testAlignmentDpFormulaComputeScoreUpperBandLinear()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigDefault> TDPValue;
    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 6);
    setScore(dpValue3, 16);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionUpperBand()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 6);

    dpValue1._score = 6;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionUpperBand()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 4);

    dpValue1._score = 4;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionUpperBand()),(+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL));
    SEQAN_ASSERT_EQ(dpValue._score, 2);

    dpValue2._score = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionUpperBand()),+TraceBitMask::HORIZONTAL);
    SEQAN_ASSERT_EQ(dpValue._score, 6);
}

void testAlignmentDpFormulaComputeScoreUpperBandLinearForbidden()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue;
    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1(true);
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionUpperBand()),+TraceBitMask::HORIZONTAL);
    SEQAN_ASSERT_EQ(dpValue._score, 2);

    dpValue2._forbidden = true;
    dpValue1._score = 2;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionUpperBand()),+TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
}

void testAlignmentDpFormulaComputeScoreUpperBandAffine()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigDefault> TDPValue;
    Score<int, Simple> score(2,-2,-4, -6);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionUpperBand()),(+TraceBitMask::DIAGONAL | TraceBitMask::HORIZONTAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);

    dpValue1._score = 6;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionUpperBand()),(+TraceBitMask::DIAGONAL | TraceBitMask::HORIZONTAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);

    dpValue2._score = 10;
    dpValue2._scoreHorizontal = 6;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionUpperBand()),(+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);

    dpValue1._score = 2;
    dpValue2._score = 6;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionUpperBand()),(+TraceBitMask::HORIZONTAL));
    SEQAN_ASSERT_EQ(dpValue._score, 2);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 2);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);

    dpValue1._score = 4;
    dpValue2._score = 6;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionUpperBand()),(+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL));
    SEQAN_ASSERT_EQ(dpValue._score, 2);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 2);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
}

void testAlignmentDpFormulaComputeScoreUpperBandAffineForbidden()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue;
    Score<int, Simple> score(2,-2,-4,-6);

    TDPValue dpValue;
    TDPValue dpValue1(true);
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionUpperBand()),+TraceBitMask::HORIZONTAL_OPEN);
    SEQAN_ASSERT_EQ(dpValue._score, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);

    dpValue2._score = 10;
    dpValue2._scoreHorizontal = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionUpperBand()),+TraceBitMask::HORIZONTAL);
    SEQAN_ASSERT_EQ(dpValue._score, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);

    dpValue1._score = 0;
    dpValue1._scoreHorizontal = 0;
    dpValue1._scoreVertical = 0;
    dpValue2._forbidden = true;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionUpperBand()),+TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);
}

void testAlignmentDpFormulaComputeScoreLowerBandLinear()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigDefault> TDPValue;
    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 16);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionLowerBand()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 6);

    dpValue1 = dpValue;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionLowerBand()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 4);

    dpValue1 = dpValue;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionLowerBand()),(+TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL));
    SEQAN_ASSERT_EQ(dpValue._score, 2);

    dpValue3._score = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionLowerBand()),+TraceBitMask::VERTICAL);
    SEQAN_ASSERT_EQ(dpValue._score, 6);
}

void testAlignmentDpFormulaComputeScoreLowerBandLinearForbidden()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue;
    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1(true);
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionLowerBand()),+TraceBitMask::VERTICAL);
    SEQAN_ASSERT_EQ(dpValue._score, 2);

    dpValue3._forbidden = true;
    dpValue1._score = dpValue._score;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionLowerBand()),+TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
}

void testAlignmentDpFormulaComputeScoreLowerBandAffine()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigDefault> TDPValue;
    Score<int, Simple> score(2,-2,-4, -6);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionLowerBand()),(+TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);

    dpValue1 = dpValue;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionLowerBand()),(+TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);

    dpValue1 = dpValue;
    dpValue3._score = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionLowerBand()),(+TraceBitMask::VERTICAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);

    dpValue1 = dpValue;
    dpValue3._score = 10;
    dpValue3._scoreVertical = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionLowerBand()),(+TraceBitMask::VERTICAL));
    SEQAN_ASSERT_EQ(dpValue._score, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 6);

    dpValue1 = dpValue;
    dpValue3._score = 8;
    dpValue3._scoreVertical = 8;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionLowerBand()),(+TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);
}

void testAlignmentDpFormulaComputeScoreLowerBandAffineForbidden()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue;
    Score<int, Simple> score(2,-2,-4, -6);

    TDPValue dpValue;
    TDPValue dpValue1(true);
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionLowerBand()),+TraceBitMask::VERTICAL_OPEN);
    SEQAN_ASSERT_EQ(dpValue._score, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 0);

    dpValue3._score = 10;
    dpValue3._scoreVertical = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionLowerBand()),+TraceBitMask::VERTICAL);
    SEQAN_ASSERT_EQ(dpValue._score, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 6);

    dpValue1._score = 0;
    dpValue3._forbidden = true;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionLowerBand()),+TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 0);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, MinValue<int>::VALUE);
}

void testAlignmentDpFormulaComputeScoreAllLinear()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigDefault> TDPValue;

    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionAll()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 6);

    dpValue1 = dpValue;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionAll()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 4);

    dpValue1 = dpValue;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionAll()),(+TraceBitMask::VERTICAL | +TraceBitMask::HORIZONTAL | +TraceBitMask::DIAGONAL));
    SEQAN_ASSERT_EQ(dpValue._score, 2);

    dpValue1 = dpValue;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'A', score, DPDirectionAll()),(+TraceBitMask::VERTICAL | +TraceBitMask::HORIZONTAL));
    SEQAN_ASSERT_EQ(dpValue._score, 2);

    dpValue1 = dpValue;
    dpValue3._score = 8;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'A', score, DPDirectionAll()),(+TraceBitMask::VERTICAL));
    SEQAN_ASSERT_EQ(dpValue._score, 4);

    dpValue1 = dpValue;
    dpValue2._score = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'A', score, DPDirectionAll()),(+TraceBitMask::HORIZONTAL));
    SEQAN_ASSERT_EQ(dpValue._score, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'C', score, DPDirectionAll()),(+TraceBitMask::HORIZONTAL | +TraceBitMask::DIAGONAL));
    SEQAN_ASSERT_EQ(dpValue._score, 6);

    dpValue1 = dpValue;
    dpValue3._score = 12;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'C', score, DPDirectionAll()),(+TraceBitMask::VERTICAL | +TraceBitMask::DIAGONAL));
    SEQAN_ASSERT_EQ(dpValue._score, 8);
}

void testAlignmentDpFormulaComputeScoreAllLinearForbidden()
{
    using namespace seqan;

    typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue;

    Score<int, Simple> score(2,-2,-4);

    TDPValue dpValue;
    TDPValue dpValue1(true);
    TDPValue dpValue2(true);
    TDPValue dpValue3(true);

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionAll()),+TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);

    dpValue1._forbidden = false;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionAll()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 6);

    dpValue1._forbidden = true;
    dpValue2._forbidden = false;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionAll()),+TraceBitMask::HORIZONTAL);
    SEQAN_ASSERT_EQ(dpValue._score, 2);

    dpValue2._forbidden = true;
    dpValue3._forbidden = false;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionAll()),+TraceBitMask::VERTICAL);
    SEQAN_ASSERT_EQ(dpValue._score, 2);
}

void testAlignmentDpFormulaComputeScoreAllAffine()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigDefault> TDPValue;
    Score<int, Simple> score(2,-2,-4, -6);

    TDPValue dpValue;
    TDPValue dpValue1;
    TDPValue dpValue2;
    TDPValue dpValue3;

    // Diagonal

    dpValue1._score = 6;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionAll()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 8);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, -4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, -4);

    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionAll()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, -4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, -4);

    // Horizontal gapOpen

    dpValue1._score = 0;
    dpValue2._scoreHorizontal = 4;
    dpValue2._score = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionAll()),+TraceBitMask::HORIZONTAL_OPEN);
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, -4);

    //Diagonal Horizontal gapOpen
    dpValue1._score = 2;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'C', score, DPDirectionAll()),(+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, -4);

    // Horizontal gapExtend
    dpValue2._scoreHorizontal = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'C', score, DPDirectionAll()),(+TraceBitMask::HORIZONTAL));
    SEQAN_ASSERT_EQ(dpValue._score, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, -4);

    //Diagonal Horizontal gapExtend
    dpValue1._score = 4;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'C', score, DPDirectionAll()),(+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL));
    SEQAN_ASSERT_EQ(dpValue._score, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, -4);

    // Vertical gapOpen
    dpValue2._score = 0;
    dpValue2._scoreHorizontal = 0;
    dpValue3._scoreVertical = 4;
    dpValue3._score = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'A', score, DPDirectionAll()),(+TraceBitMask::VERTICAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, -4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);

    // Vertical gapExtend
    dpValue3._scoreVertical = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'A', score, DPDirectionAll()),(+TraceBitMask::VERTICAL));
    SEQAN_ASSERT_EQ(dpValue._score, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, -4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 6);

    //Diagonal Vertical gapOpen
    dpValue1._score = 2;
    dpValue3._scoreVertical = 4;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'C', score, DPDirectionAll()),(+TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, -4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);

    //Diagonal Vertical gapExtend
    dpValue1._score = 4;
    dpValue3._scoreVertical = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'C', score, DPDirectionAll()),(+TraceBitMask::VERTICAL | +TraceBitMask::DIAGONAL));
    SEQAN_ASSERT_EQ(dpValue._score, 6);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, -4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 6);

    //Horizontal gapOpen Vertical gapOpen
    dpValue2._score = 10;
    dpValue2._scoreHorizontal = 4;
    dpValue3._scoreVertical = 4;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'G', score, DPDirectionAll()),(+TraceBitMask::VERTICAL_OPEN | +TraceBitMask::HORIZONTAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);

    //Horizontal gapExtend Vertical gapOpen
    dpValue2._scoreHorizontal = 8;
    dpValue2._score = 8;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'G', score, DPDirectionAll()),(+TraceBitMask::VERTICAL_OPEN | +TraceBitMask::HORIZONTAL));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);

    //Horizontal gapOpen Vertical gapExtend
    dpValue2._scoreHorizontal = 4;
    dpValue2._score = 10;
    dpValue3._scoreVertical = 8;
    dpValue3._score = 8;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'G', 'A', score, DPDirectionAll()),(+TraceBitMask::VERTICAL | +TraceBitMask::HORIZONTAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);

    //Horizontal gapExtend Vertical gapExtend
    dpValue2._scoreHorizontal = 8;
    dpValue2._score = 8;
    dpValue3._scoreVertical = 8;
    dpValue3._score = 8;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'G', 'A', score, DPDirectionAll()),(+TraceBitMask::VERTICAL | +TraceBitMask::HORIZONTAL));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);

    // Diagonal Horizontal gapOpen Vertical gapOpen
    dpValue1._score = 2;
    dpValue2._scoreHorizontal = 4;
    dpValue2._score = 10;
    dpValue3._scoreVertical = 4;
    dpValue3._score = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'G', 'G', score, DPDirectionAll()),(+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL_OPEN | +TraceBitMask::VERTICAL_OPEN));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);

    // Diagonal Horizontal gapExtend Vertical gapOpen
    dpValue2._scoreHorizontal = 8;
    dpValue2._score = 8;
    dpValue3._scoreVertical = 4;
    dpValue3._score = 10;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'G', 'G', score, DPDirectionAll()),(+TraceBitMask::HORIZONTAL | +TraceBitMask::VERTICAL_OPEN | +TraceBitMask::DIAGONAL));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);

    // Diagonal Horizontal gapOpen Vertical gapExtend
    dpValue2._scoreHorizontal = 4;
    dpValue2._score = 10;
    dpValue3._scoreVertical = 8;
    dpValue3._score = 8;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'G', 'G', score, DPDirectionAll()),(+TraceBitMask::HORIZONTAL_OPEN | +TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);

    // Diagonal Horizontal gapExtend Vertical gapExtend
    dpValue2._scoreHorizontal = 8;
    dpValue2._score = 8;
    dpValue3._scoreVertical = 8;
    dpValue3._score = 8;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'G', 'G', score, DPDirectionAll()),(+TraceBitMask::HORIZONTAL | +TraceBitMask::DIAGONAL | +TraceBitMask::VERTICAL));
    SEQAN_ASSERT_EQ(dpValue._score, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreHorizontal, 4);
    SEQAN_ASSERT_EQ(dpValue._scoreVertical, 4);
}

void testAlignmentDpFormulaComputeScoreAllAffineForbidden()
{
    using namespace seqan;

    typedef DPValue<int, AffineGaps, DPValueConfigForbidden> TDPValue;

    Score<int, Simple> score(2,-2,-4, -6);

    TDPValue dpValue;
    TDPValue dpValue1(true);
    TDPValue dpValue2(true);
    TDPValue dpValue3(true);

    setScore(dpValue, 6);
    setScore(dpValue1, 4);
    setScore(dpValue2, 6);
    setScore(dpValue3, 6);


    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionAll()),+TraceBitMask::FORBIDDEN);
    SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);

    dpValue1._forbidden = false;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionAll()),+TraceBitMask::DIAGONAL);
    SEQAN_ASSERT_EQ(dpValue._score, 6);

    dpValue1._forbidden = true;
    dpValue2._forbidden = false;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionAll()),+TraceBitMask::HORIZONTAL_OPEN);
    SEQAN_ASSERT_EQ(dpValue._score, 0);

    dpValue2._forbidden = true;
    dpValue3._forbidden = false;
    SEQAN_ASSERT_EQ(computeScore(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', score, DPDirectionAll()),+TraceBitMask::VERTICAL_OPEN);
    SEQAN_ASSERT_EQ(dpValue._score, 0);
}

void testAlignmentDpFormulaGlobal()
{
    using namespace seqan;

    { // forbidden disabled
        typedef DPValue<int, LinearGaps, DPValueConfigDefault> TDPValue;
        typedef Score<int, Simple> TScoreScheme;
        typedef DPFormula<TScoreScheme, Global<> > TDPFormula;

        TDPFormula dpFormula(TScoreScheme(2,-2,-4));

        TDPValue dpValue;
        TDPValue dpValue1;
        TDPValue dpValue2;
        TDPValue dpValue3;

        setScore(dpValue, 6);
        setScore(dpValue1, 4);
        setScore(dpValue2, 6);
        setScore(dpValue3, 6);

        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),+TraceBitMask::DIAGONAL);
        SEQAN_ASSERT_EQ(dpValue._score, 6);
    }


    { // forbidden enabled
        typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue;
        typedef Score<int, Simple> TScoreScheme;
        typedef DPFormula<TScoreScheme, Global<> > TDPFormula;

        TDPFormula dpFormula(TScoreScheme(2,-2,-4));

        TDPValue dpValue;
        TDPValue dpValue1;
        TDPValue dpValue2;
        TDPValue dpValue3;

        setScore(dpValue, 6);
        setScore(dpValue1, 4);
        setScore(dpValue2, 6);
        setScore(dpValue3, 6);

        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),+TraceBitMask::DIAGONAL);
        SEQAN_ASSERT_EQ(dpValue._score, 6);

        dpValue1._forbidden = true;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),(+TraceBitMask::HORIZONTAL | +TraceBitMask::VERTICAL));
        SEQAN_ASSERT_EQ(dpValue._score, 2);

        dpValue2._forbidden = true;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),(+TraceBitMask::VERTICAL));
        SEQAN_ASSERT_EQ(dpValue._score, 2);

        dpValue3._forbidden = true;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),(+TraceBitMask::FORBIDDEN));
        SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
    }
}

void testAlignmentDpFormulaLocal()
{
    using namespace seqan;

    { // forbidden disabled
        typedef DPValue<int, LinearGaps, DPValueConfigDefault> TDPValue;
        typedef Score<int, Simple> TScoreScheme;
        typedef DPFormula<TScoreScheme, Local<> > TDPFormula;

        TDPFormula dpFormula(TScoreScheme(2,-2,-4));

        TDPValue dpValue;
        TDPValue dpValue1;
        TDPValue dpValue2;
        TDPValue dpValue3;

        setScore(dpValue, 6);
        setScore(dpValue1, 4);
        setScore(dpValue2, 6);
        setScore(dpValue3, 6);

        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),+TraceBitMask::DIAGONAL);
        SEQAN_ASSERT_EQ(dpValue._score, 6);

        dpValue1._score = 2;
        dpValue2._score = 4;
        dpValue3._score = 4;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'A', DPDirectionAll()),+TraceBitMask::NONE);
        SEQAN_ASSERT_EQ(dpValue._score, 0);

        dpValue1._score = 3;
        dpValue2._score = 5;
        dpValue3._score = 5;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'A', DPDirectionAll()),
                        (+TraceBitMask::DIAGONAL | +TraceBitMask::HORIZONTAL | + TraceBitMask::VERTICAL));
        SEQAN_ASSERT_EQ(dpValue._score, 1);

        dpValue1._score = 0;
        dpValue2._score = 5;
        dpValue3._score = 2;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'A', DPDirectionAll()),
                        (+TraceBitMask::HORIZONTAL));
        SEQAN_ASSERT_EQ(dpValue._score, 1);

        dpValue1._score = -10;
        dpValue2._score = -2;
        dpValue3._score = -2;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),
                        (+TraceBitMask::NONE));
        SEQAN_ASSERT_EQ(dpValue._score, 0);
    }


    { // forbidden enabled
        typedef DPValue<int, LinearGaps, DPValueConfigForbidden> TDPValue;
        typedef Score<int, Simple> TScoreScheme;
        typedef DPFormula<TScoreScheme, Local<> > TDPFormula;

        TDPFormula dpFormula(TScoreScheme(2,-2,-4));

        TDPValue dpValue;
        TDPValue dpValue1;
        TDPValue dpValue2;
        TDPValue dpValue3;

        setScore(dpValue, 6);
        setScore(dpValue1, 4);
        setScore(dpValue2, 6);
        setScore(dpValue3, 6);

        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),+TraceBitMask::DIAGONAL);
        SEQAN_ASSERT_EQ(dpValue._score, 6);

        dpValue1._score = 2;
        dpValue2._score = 4;
        dpValue3._score = 4;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'C', 'A', DPDirectionAll()),+TraceBitMask::NONE);
        SEQAN_ASSERT_EQ(dpValue._score, 0);

        dpValue1._forbidden = true;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),(+TraceBitMask::NONE));
        SEQAN_ASSERT_EQ(dpValue._score, 0);

        dpValue2._forbidden = true;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),(+TraceBitMask::NONE));
        SEQAN_ASSERT_EQ(dpValue._score, 0);

        dpValue3._forbidden = true;
        SEQAN_ASSERT_EQ(dpFormula(dpValue, dpValue1, dpValue2, dpValue3, 'A', 'A', DPDirectionAll()),(+TraceBitMask::FORBIDDEN));
        SEQAN_ASSERT_EQ(dpValue._score, MinValue<int>::VALUE);
    }
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_linear_diagonal)
{
    testAlignmentDpFormulaComputeScoreDiagonalLinear();
    testAlignmentDpFormulaComputeScoreDiagonalLinearForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_linear_vertical)
{
    testAlignmentDpFormulaComputeScoreVerticalLinear();
    testAlignmentDpFormulaComputeScoreVerticalLinearForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_linear_horizontal)
{
    testAlignmentDpFormulaComputeScoreHorizontalLinear();
    testAlignmentDpFormulaComputeScoreHorizontalLinearForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_linear_upper_band)
{
    testAlignmentDpFormulaComputeScoreUpperBandLinear();
    testAlignmentDpFormulaComputeScoreUpperBandLinearForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_linear_lower_band)
{
    testAlignmentDpFormulaComputeScoreLowerBandLinear();
    testAlignmentDpFormulaComputeScoreLowerBandLinearForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_linear_all)
{
    testAlignmentDpFormulaComputeScoreAllLinear();
    testAlignmentDpFormulaComputeScoreAllLinearForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_linear_none)
{
    testAlignmentDpFormulaComputeScoreNoneLinear();
    testAlignmentDpFormulaComputeScoreNoneLinearForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_affine_diagonal)
{
    testAlignmentDpFormulaComputeScoreDiagonalAffine();
    testAlignmentDpFormulaComputeScoreDiagonalAffineForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_affine_vertical)
{
    testAlignmentDpFormulaComputeScoreVerticalAffine();
    testAlignmentDpFormulaComputeScoreVerticalAffineForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_affine_horizontal)
{
    testAlignmentDpFormulaComputeScoreHorizontalAffine();
    testAlignmentDpFormulaComputeScoreHorizontalAffineForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_affine_upper_band)
{
    testAlignmentDpFormulaComputeScoreUpperBandAffine();
    testAlignmentDpFormulaComputeScoreUpperBandAffineForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_affine_lower_band)
{
    testAlignmentDpFormulaComputeScoreLowerBandAffine();
    testAlignmentDpFormulaComputeScoreLowerBandAffineForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_affine_all)
{
    testAlignmentDpFormulaComputeScoreAllAffine();
    testAlignmentDpFormulaComputeScoreAllAffineForbidden();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_affine_none)
{
    testAlignmentDpFormulaComputeScoreNoneAffine();
    testAlignmentDpFormulaComputeScoreNoneAffineForbidden();
}

    // dp formula object
SEQAN_DEFINE_TEST(test_alignment_dp_formula_global)
{
    testAlignmentDpFormulaGlobal();
}

SEQAN_DEFINE_TEST(test_alignment_dp_formula_local)
{
    testAlignmentDpFormulaLocal();
}

#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_FORMULA_H_
