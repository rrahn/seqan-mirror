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

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_MATRIX_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_MATRIX_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/align/alignment_base.h>
#include <seqan/align/alignment_dp_band.h>
#include <seqan/align/alignment_dp_value.h>
#include <seqan/align/alignment_dp_matrix.h>

using namespace seqan;

template <typename TDPValue>
void testAlignmentDPMatrixConstructor()
{
    DPMatrix<TDPValue> dpMatrix;

    SEQAN_ASSERT_EQ(length(dpMatrix._dataTable), 0u);

    CharString str0 = "Hello World";
    CharString str1 = "Hello Seqan";

    DPMatrix<TDPValue> dpMatrix2(str0, str1);
    SEQAN_ASSERT_EQ(length(dpMatrix2._dataTable), 12u);

    DPMatrix<TDPValue> dpMatrix3(str0, str1, Band<BandSwitchedOff>());
    SEQAN_ASSERT_EQ(length(dpMatrix3._dataTable), 12u);

    DPMatrix<TDPValue> dpMatrix4(str0, str1, Band<BandSwitchedOn<> >(-2,2));
    SEQAN_ASSERT_EQ(length(dpMatrix4._dataTable), 5u);

    DPMatrix<TDPValue> dpMatrix5(str0, str1, Band<BandSwitchedOn<> >(2,2));
    SEQAN_ASSERT_EQ(length(dpMatrix5._dataTable), 1u);
}

template <typename TDPValue>
void testAlignmentDPMatrixInitDpMatrix()
{
    DPMatrix<TDPValue> dpMatrix;

    SEQAN_ASSERT_EQ(length(dpMatrix._dataTable), 0u);

    CharString str0 = "Hello World";
    CharString str1 = "Hello Seqan";

    _initDPMatrix(dpMatrix, str0, str1, Band<BandSwitchedOff>());

    SEQAN_ASSERT_EQ(length(dpMatrix._dataTable), 12u);

    _initDPMatrix(dpMatrix, str0, str1, Band<BandSwitchedOn<> >(-2,2));
    SEQAN_ASSERT_EQ(length(dpMatrix._dataTable), 5u);

    _initDPMatrix(dpMatrix, str0, str1, Band<BandSwitchedOn<> >(2,2));
    SEQAN_ASSERT_EQ(length(dpMatrix._dataTable), 1u);
}

template <typename TDPValue>
void testAlignmentDPMatrixDpMatrixValue()
{
    typedef DPMatrix<TDPValue> TDPMatrix;
    typedef DPMatrix<TDPValue> const TConstDPMatrix;
    typedef typename IsSameType<TDPValue, typename Value<TDPMatrix>::Type >::Type TResult1;
    typedef typename IsSameType<TDPValue const, typename Value<TConstDPMatrix>::Type >::Type TResult2;
    SEQAN_ASSERT_EQ(+TResult1::VALUE, true);
    SEQAN_ASSERT_EQ(+TResult2::VALUE, true);
}

template <typename TDPValue>
void testAlignmentDPMatrixDpMatrixReference()
{
    typedef DPMatrix<TDPValue> TDPMatrix;
    typedef DPMatrix<TDPValue> const TConstDPMatrix;
    typedef typename IsSameType<TDPValue &, typename Reference<TDPMatrix>::Type >::Type TResult1;
    typedef typename IsSameType<TDPValue const & , typename Reference<TConstDPMatrix>::Type >::Type TResult2;
    SEQAN_ASSERT_EQ(+TResult1::VALUE, true);
    SEQAN_ASSERT_EQ(+TResult2::VALUE, true);
}

template <typename TDPValue>
void testAlignmentDPMatrixDpMatrixSize()
{
    typedef DPMatrix<TDPValue> TDPMatrix;
    typedef typename TDPMatrix::TDataTable_ TDataTable_;
    typedef typename TDPMatrix::TDataTable_ const TConstDataTable_;
    typedef DPMatrix<TDPValue> const TConstDPMatrix;
    typedef typename IsSameType<typename Size<TDataTable_>::Type, typename Size<TDPMatrix>::Type >::Type TResult1;
    typedef typename IsSameType<typename Size<TConstDataTable_>::Type, typename Size<TConstDPMatrix>::Type >::Type TResult2;
    SEQAN_ASSERT_EQ(+TResult1::VALUE, true);
    SEQAN_ASSERT_EQ(+TResult2::VALUE, true);
}

template <typename TDPValue>
void testAlignmentDPMatrixDpMatrixPosition()
{
    typedef DPMatrix<TDPValue> TDPMatrix;
    typedef typename TDPMatrix::TDataTable_ TDataTable_;
    typedef typename TDPMatrix::TDataTable_ const TConstDataTable_;
    typedef DPMatrix<TDPValue> const TConstDPMatrix;

    typedef typename IsSameType<typename Position<TDataTable_>::Type, typename Position<TDPMatrix>::Type >::Type TResult1;
    typedef typename IsSameType<typename Position<TConstDataTable_>::Type, typename Position<TConstDPMatrix>::Type >::Type TResult2;

    SEQAN_ASSERT_EQ(+TResult1::VALUE, true);
    SEQAN_ASSERT_EQ(+TResult2::VALUE, true);
}

template <typename TDPValue>
void testAlignmentDPMatrixDpMatrixIterator()
{
    typedef DPMatrix<TDPValue> TDPMatrix;
    typedef DPMatrix<TDPValue> const TConstDPMatrix;
    typedef typename TDPMatrix::TDataTable_ TDataTable_;
    typedef typename TDPMatrix::TDataTable_ const TConstDataTable_;
    typedef typename Iterator<TDPMatrix, Standard>::Type TIterator;
    typedef typename Iterator<TConstDPMatrix, Standard>::Type TConstIterator;

    typedef typename IsSameType<typename Iterator<TDataTable_>::Type, TIterator >::Type TResult1;
    typedef typename IsSameType<typename Iterator<TConstDataTable_>::Type, TConstIterator >::Type TResult2;

    SEQAN_ASSERT_EQ(+TResult1::VALUE, true);
    SEQAN_ASSERT_EQ(+TResult2::VALUE, true);
}

template <typename TDPValue>
void testAlignmentDPMatrixClear()
{
    CharString str0 = "Hello Seqan";
    CharString str1 = "Hello Team";

    DPMatrix<TDPValue> dpMatrix(str0, str1);
    SEQAN_ASSERT_EQ(length(dpMatrix._dataTable), length(str1)+1);
    clear(dpMatrix);
    SEQAN_ASSERT_EQ(length(dpMatrix._dataTable), 0u);

    DPMatrix<TDPValue> dpMatrix1(str0, str1, Band<BandSwitchedOn<> >(-2,2));
    SEQAN_ASSERT_EQ(length(dpMatrix1._dataTable), 5u);
    clear(dpMatrix1);
    SEQAN_ASSERT_EQ(length(dpMatrix1._dataTable), 0u);

    DPMatrix<TDPValue> dpMatrix2(str0, str1, Band<BandSwitchedOn<> >(-2,2));
    SEQAN_ASSERT_EQ(length(dpMatrix2._dataTable), 5u);
    clear(dpMatrix2);
    SEQAN_ASSERT_EQ(length(dpMatrix2._dataTable), 0u);
}

template <typename TDPValue>
void testAlignmentDPMatrixLength()
{
    CharString str0 = "Hello Seqan";
    CharString str1 = "Hello Team";

    DPMatrix<TDPValue> dpMatrix(str0, str1);
    SEQAN_ASSERT_EQ(length(dpMatrix._dataTable), length(dpMatrix) + 1);

    DPMatrix<TDPValue> dpMatrix1(str0, str1, Band<BandSwitchedOn<> >(-2,2));
    SEQAN_ASSERT_EQ(length(dpMatrix1._dataTable), 3u);
    clear(dpMatrix);
    SEQAN_ASSERT_EQ(length(dpMatrix1._dataTable), length(dpMatrix1));

    DPMatrix<TDPValue> dpMatrix2(str0, str1, Band<BandSwitchedOn<> >(-2,2));
    SEQAN_ASSERT_EQ(length(dpMatrix2._dataTable), 3u);
    clear(dpMatrix);
    SEQAN_ASSERT_EQ(length(dpMatrix2._dataTable), length(dpMatrix2));
}


template <typename TDPValue>
void testAlignmentDPMatrixBegin()
{
    typedef DPMatrix<TDPValue> TDPMatrix;
    typedef typename Iterator<TDPMatrix, Standard>::Type TIterator;

    CharString str0 = "Hello Seqan";
    CharString str1 = "Hello Team";

    TDPMatrix dpMatrix(str0, str1);
    setScore(dpMatrix._dataTable[0], 10);
    setScore(dpMatrix._dataTable[1], -20);

    TIterator dpBegin = begin(dpMatrix);
    SEQAN_ASSERT_EQ(getScore(value(dpBegin)), 10);
    SEQAN_ASSERT_EQ(getScore(value(++dpBegin)), -20);

    TDPMatrix dpMatrix1(str0, str1, Band<BandSwitchedOn<> >(-2,3));
    setScore(dpMatrix1._dataTable[0], 10);
    setScore(dpMatrix1._dataTable[1], -20);

    dpBegin = begin(dpMatrix1);
    SEQAN_ASSERT_EQ(getScore(value(dpBegin)), 10);
    SEQAN_ASSERT_EQ(getScore(value(++dpBegin)), -20);

    TDPMatrix dpMatrix2(str0, str1, Band<BandSwitchedOn<> >(-2,3));
    setScore(dpMatrix2._dataTable[0], 10);
    setScore(dpMatrix2._dataTable[1], -20);

    dpBegin = begin(dpMatrix2);
    SEQAN_ASSERT_EQ(getScore(value(dpBegin)), 10);
    SEQAN_ASSERT_EQ(getScore(value(++dpBegin)), -20);
}


template <typename TDPValue>
void testAlignmentDPMatrixEnd()
{
    typedef DPMatrix<TDPValue> TDPMatrix;
    typedef typename Iterator<TDPMatrix, Standard>::Type TIterator;

    CharString str0 = "Hello Seqan";
    CharString str1 = "Hello Team";

    TDPMatrix dpMatrix(str0, str1);
    setScore(dpMatrix._dataTable[length(str1)], 10);
    setScore(dpMatrix._dataTable[length(str1)-1], -20);

    TIterator dpEnd = end(dpMatrix);
    SEQAN_ASSERT_EQ(getScore(value(--dpEnd)), 10);
    SEQAN_ASSERT_EQ(getScore(value(--dpEnd)), -20);

    TDPMatrix dpMatrix1(str0, str1, Band<BandSwitchedOn<> >(-2,3));
    setScore(dpMatrix1._dataTable[5], 10);
    setScore(dpMatrix1._dataTable[4], -20);

    dpEnd = end(dpMatrix1);
    SEQAN_ASSERT_EQ(getScore(value(--dpEnd)), 10);
    SEQAN_ASSERT_EQ(getScore(value(--dpEnd)), -20);

    TDPMatrix dpMatrix2(str0, str1, Band<BandSwitchedOn<> >(-2,3));
    setScore(dpMatrix2._dataTable[5], 10);
    setScore(dpMatrix2._dataTable[4], -20);

    dpEnd = end(dpMatrix2);
    SEQAN_ASSERT_EQ(getScore(value(--dpEnd)), 10);
    SEQAN_ASSERT_EQ(getScore(value(--dpEnd)), -20);
}

template<typename TDPValue>
void testAlignmentDPMatrixGetColumnSize()
{

    CharString str0 = "Hello Seqan";
    CharString str1 = "Hello Team";


    SEQAN_ASSERT_EQ(getColumnSize(str0, str1, Band<BandSwitchedOff>()), length(str1) + 1);
    SEQAN_ASSERT_EQ(getColumnSize(str0, str1, Band<BandSwitchedOn<> >(-2,2)), 5u);
    SEQAN_ASSERT_EQ(getColumnSize(str0, str1, Band<BandSwitchedOn<> >(-2,3)), 6u);
    SEQAN_ASSERT_EQ(getColumnSize(str0, str1, Band<BandSwitchedOn<> >(-4,3)), 8u);
}

template<typename TDPValue>
void testAlignmentDPMatrixGetRowSize()
{

    CharString str0 = "Hello Seqan";
    CharString str1 = "Hello Team";


    SEQAN_ASSERT_EQ(getRowSize(str0, str1, Band<BandSwitchedOff>()), length(str0) + 1);
    SEQAN_ASSERT_EQ(getRowSize(str0, str1, Band<BandSwitchedOn<> >(-2,3)), 6u);
    SEQAN_ASSERT_EQ(getRowSize(str0, str1, Band<BandSwitchedOn<> >(-4,3)), 8u);
}

template <typename TDPValue>
void testAlignmentDPMatrixGetIntitialColumnSize()
{
    CharString str0 = "Hello Seqan";
    CharString str1 = "Hello Team";

    SEQAN_ASSERT_EQ(getInitialColumnSize(str0, str1, Band<BandSwitchedOff>()), length(str1) + 1);
    SEQAN_ASSERT_EQ(getInitialColumnSize(str0, str1, Band<BandSwitchedOn<> >(-2,2)), 3u);
    SEQAN_ASSERT_EQ(getInitialColumnSize(str0, str1, Band<BandSwitchedOn<> >(-5,2)), 6u);
    SEQAN_ASSERT_EQ(getInitialColumnSize(str0, str1, Band<BandSwitchedOn<> >(-76,-5)), 6u);
    SEQAN_ASSERT_EQ(getInitialColumnSize(str0, str1, Band<BandSwitchedOn<> >(-76,-28)), 1u);
    SEQAN_ASSERT_EQ(getInitialColumnSize(str0, str1, Band<BandSwitchedOn<> >(28, 76)), 1u);
}


template <typename TDPValue>
void testAlignmentDPMatrixGetIntitialRowSize()
{
    CharString str0 = "Hello Seqan";
    CharString str1 = "Hello Team";

    SEQAN_ASSERT_EQ(getInitialRowSize(str0, str1, Band<BandSwitchedOff>()), length(str0) + 1);
    SEQAN_ASSERT_EQ(getInitialRowSize(str0, str1, Band<BandSwitchedOn<> >(-2,2)), 3u);
    SEQAN_ASSERT_EQ(getInitialRowSize(str0, str1, Band<BandSwitchedOn<> >(-2,5)), 6u);
    SEQAN_ASSERT_EQ(getInitialRowSize(str0, str1, Band<BandSwitchedOn<> >(28, 76)), 12u);
    SEQAN_ASSERT_EQ(getInitialRowSize(str0, str1, Band<BandSwitchedOn<> >(-76,-28)), 1u);
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_constructor)
{
    testAlignmentDPMatrixConstructor<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixConstructor<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_init_dp_matrix)
{
    testAlignmentDPMatrixInitDpMatrix<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixInitDpMatrix<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_value)
{
    testAlignmentDPMatrixDpMatrixValue<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixDpMatrixValue<DPValue<int, AffineGaps> >();
}


SEQAN_DEFINE_TEST(test_alignment_dp_matrix_reference)
{
    testAlignmentDPMatrixDpMatrixReference<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixDpMatrixReference<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_size)
{
    testAlignmentDPMatrixDpMatrixSize<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixDpMatrixSize<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_position)
{
    testAlignmentDPMatrixDpMatrixPosition<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixDpMatrixPosition<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_iterator)
{
    testAlignmentDPMatrixDpMatrixIterator<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixDpMatrixIterator<DPValue<int, AffineGaps> >();
}


SEQAN_DEFINE_TEST(test_alignment_dp_matrix_clear)
{
    testAlignmentDPMatrixClear<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixClear<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_length)
{
    testAlignmentDPMatrixClear<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixClear<DPValue<int, AffineGaps> >();
}


SEQAN_DEFINE_TEST(test_alignment_dp_matrix_begin)
{
    testAlignmentDPMatrixBegin<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixBegin<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_end)
{
    testAlignmentDPMatrixEnd<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixEnd<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_get_column_size)
{
    testAlignmentDPMatrixGetColumnSize<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixGetColumnSize<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_get_row_size)
{
    testAlignmentDPMatrixGetRowSize<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixGetRowSize<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_get_initial_column_size)
{
    testAlignmentDPMatrixGetIntitialColumnSize<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixGetIntitialColumnSize<DPValue<int, AffineGaps> >();
}

SEQAN_DEFINE_TEST(test_alignment_dp_matrix_get_initial_row_size)
{
    testAlignmentDPMatrixGetIntitialRowSize<DPValue<int, LinearGaps> >();
    testAlignmentDPMatrixGetIntitialRowSize<DPValue<int, AffineGaps> >();
}

#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_MATRIX_H_
