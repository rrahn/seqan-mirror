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

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_BAND_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_BAND_H_

#include <seqan/sequence.h>

#include <seqan/align/alignment_base.h>
#include <seqan/align/alignment_dp_band.h>

using namespace seqan;

void testAlignmentDpBandConstructor()
{
    { // band off
        Band<BandSwitchedOff> band;

        SEQAN_ASSERT_EQ(band._upperDiagonal, 0);
        SEQAN_ASSERT_EQ(band._lowerDiagonal, 0);

    }

    { // default band
        Band<BandSwitchedOn<> > band;

        SEQAN_ASSERT_EQ(band._upperDiagonal, 0);
        SEQAN_ASSERT_EQ(band._lowerDiagonal, 0);

        band._upperDiagonal = 10;
        band._lowerDiagonal = 12;

        SEQAN_ASSERT_EQ(band._upperDiagonal, 10);
        SEQAN_ASSERT_EQ(band._lowerDiagonal, 12);

        Band<BandSwitchedOn<> > band2a(3, 7);

        SEQAN_ASSERT_EQ(band2a._upperDiagonal, 7);
        SEQAN_ASSERT_EQ(band2a._lowerDiagonal, 3);

        Band<BandSwitchedOn<> > band2b(-7, -3);

        SEQAN_ASSERT_EQ(band2b._upperDiagonal, -3);
        SEQAN_ASSERT_EQ(band2b._lowerDiagonal, -7);

        Band<BandSwitchedOn<> > band3(band2b);

        SEQAN_ASSERT_EQ(band3._upperDiagonal, -3);
        SEQAN_ASSERT_EQ(band3._lowerDiagonal, -7);
    }

    { // wide band
        Band<BandSwitchedOn<WideBand> > band;

        SEQAN_ASSERT_EQ(band._upperDiagonal, 0);
        SEQAN_ASSERT_EQ(band._lowerDiagonal, 0);

        band._upperDiagonal = 10;
        band._lowerDiagonal = 12;

        SEQAN_ASSERT_EQ(band._upperDiagonal, 10);
        SEQAN_ASSERT_EQ(band._lowerDiagonal, 12);

        Band<BandSwitchedOn<> > band2a(3, 7);

        SEQAN_ASSERT_EQ(band2a._upperDiagonal, 7);
        SEQAN_ASSERT_EQ(band2a._lowerDiagonal, 3);

        Band<BandSwitchedOn<> > band2b(-7, -3);

        SEQAN_ASSERT_EQ(band2b._upperDiagonal, -3);
        SEQAN_ASSERT_EQ(band2b._lowerDiagonal, -7);

        Band<BandSwitchedOn<> > band3(band2b);

        SEQAN_ASSERT_EQ(band3._upperDiagonal, -3);
        SEQAN_ASSERT_EQ(band3._lowerDiagonal, -7);
    }
}


template <typename TBand>
void testAlignmentDpBandGetUpperDiagonal(TBand const & band)
{
    SEQAN_ASSERT_EQ(band._upperDiagonal, getUpperDiagonal(band));
}

template <typename TBand>
void testAlignmentDpBandGetLowerDiagonal(TBand const & band)
{
    SEQAN_ASSERT_EQ(band._lowerDiagonal, getLowerDiagonal(band));
}

template <typename TBand, typename TSize>
void testAlignmentDpBandSetUpperDiagonal(TBand & band, TSize const & horiPos)
{
    setUpperDiagonal(band, horiPos);
    SEQAN_ASSERT_EQ(band._upperDiagonal, horiPos);
}

template <typename TBand, typename TSize>
void testAlignmentDpBandSetLowerDiagonal(TBand & band, TSize const & vertiPos)
{
    setLowerDiagonal(band, vertiPos);
    SEQAN_ASSERT_EQ(band._lowerDiagonal, vertiPos);
}

template <typename TBand>
void testAlignmentDpBandGetBandSize(TBand const & band)
{
    SEQAN_ASSERT_EQ(band._upperDiagonal - band._lowerDiagonal + 1, getBandSize(band));
}

SEQAN_DEFINE_TEST(test_alignment_dp_band_constructor)
{
    testAlignmentDpBandConstructor();
}

SEQAN_DEFINE_TEST(test_alignment_dp_band_get_horizontal_pos)
{
    testAlignmentDpBandGetUpperDiagonal(Band<BandSwitchedOff>());
    testAlignmentDpBandGetUpperDiagonal(Band<BandSwitchedOn<> >(-2,7));
    testAlignmentDpBandGetUpperDiagonal(Band<BandSwitchedOn<WideBand> >(-2,7));
}

SEQAN_DEFINE_TEST(test_alignment_dp_band_get_vertical_pos)
{
    testAlignmentDpBandGetLowerDiagonal(Band<BandSwitchedOff>());
    testAlignmentDpBandGetLowerDiagonal(Band<BandSwitchedOn<> >(-2,7));
    testAlignmentDpBandGetLowerDiagonal(Band<BandSwitchedOn<WideBand> >(-2,7));
}

SEQAN_DEFINE_TEST(test_alignment_dp_band_set_horizontal_pos)
{
    Band<BandSwitchedOn<> > band(4,6);
    testAlignmentDpBandSetUpperDiagonal(band, 10);
    Band<BandSwitchedOn<> > band2;
    testAlignmentDpBandSetUpperDiagonal(band2, 4);
    Band<BandSwitchedOn<WideBand> > band3(4,6);
    testAlignmentDpBandSetUpperDiagonal(band3, 10);
    Band<BandSwitchedOn<WideBand> > band4;
    testAlignmentDpBandSetUpperDiagonal(band4, 4);
}

SEQAN_DEFINE_TEST(test_alignment_dp_band_set_vertical_pos)
{
    Band<BandSwitchedOn<> > band(4,6);
    testAlignmentDpBandSetLowerDiagonal(band, -3);
    Band<BandSwitchedOn<> > band2;
    testAlignmentDpBandSetUpperDiagonal(band2, 3);
    Band<BandSwitchedOn<WideBand> > band3(4,6);
    testAlignmentDpBandSetLowerDiagonal(band3, -3);
    Band<BandSwitchedOn<WideBand> > band4;
    testAlignmentDpBandSetUpperDiagonal(band4, 3);
}

SEQAN_DEFINE_TEST(test_alignment_dp_band_get_band_size)
{
    testAlignmentDpBandGetBandSize(Band<BandSwitchedOff>());
    testAlignmentDpBandGetBandSize(Band<BandSwitchedOn<> >(-8, -2));
    testAlignmentDpBandGetBandSize(Band<BandSwitchedOn<> >(5, 10));
    testAlignmentDpBandGetBandSize(Band<BandSwitchedOn<> >(-5, 10));
    testAlignmentDpBandGetBandSize(Band<BandSwitchedOn<WideBand> >(-8, -2));
    testAlignmentDpBandGetBandSize(Band<BandSwitchedOn<WideBand> >(5, 10));
    testAlignmentDpBandGetBandSize(Band<BandSwitchedOn<WideBand> >(-5, 10));
}

SEQAN_DEFINE_TEST(test_aligmenmt_dp_band_check_valid_band)
{
    // TODO(rmaerker): write me!
}

SEQAN_DEFINE_TEST(test_aligmenmt_dp_band_is_band_enabled)
{
    // TODO(rmaerker): write me!
}

SEQAN_DEFINE_TEST(test_aligmenmt_dp_band_is_wide_band)
{
    // TODO(rmaerker): write me!
}

#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_DP_BAND_H_
