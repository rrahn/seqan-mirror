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
// Implementation of the Band structure and its getter and setter methods.
// The band is used to restrict the dynamic programming area. In case of
// unbanded alignments the band is disabled by using the tag BandSwitschedOff.

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_BAND_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_BAND_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <>
struct Band<BandSwitchedOff>
{
    int _lowerDiagonal;
    int _upperDiagonal;

    Band() : _lowerDiagonal(0), _upperDiagonal(0){  }
};

// ----------------------------------------------------------------------------
// Class Band                                                         [Enabled]
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Band<BandSwitchedOn<TSpec> >
{
    int _lowerDiagonal;
    int _upperDiagonal;

    Band() : _lowerDiagonal(0), _upperDiagonal(0){  }

    Band(Band const & other) : _lowerDiagonal(other._lowerDiagonal), _upperDiagonal(other._upperDiagonal)
    {}

    Band(int const & lowerDiagonal, int const & upperDiagonal)
                : _lowerDiagonal(lowerDiagonal),
                     _upperDiagonal(upperDiagonal)
    {
        SEQAN_ASSERT(_checkValidBand(*this));
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TTag>
struct Position<Band<TTag> >
{
    typedef int Type;
};

template <typename TTag>
struct Position<Band<TTag> const >  : Size<Band<TTag> >{};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TTag>
struct Size<Band<TTag> >
{
    typedef int Type;
};

template <typename TTag>
struct Size<Band<TTag> const >  : Size<Band<TTag> >{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getBandSize()
// ----------------------------------------------------------------------------

template <typename TBandSpec>
inline typename Size<Band<TBandSpec> >::Type
getBandSize(Band<TBandSpec> const & band)
{
    return 1 + band._upperDiagonal - band._lowerDiagonal;  // add 1 because of 0 indexed sequences
}

// ----------------------------------------------------------------------------
// Function getHorizontalPos()
// ----------------------------------------------------------------------------

template <typename TBandSpec>
inline typename Position<Band<TBandSpec> >::Type
getUpperDiagonal(Band<TBandSpec> const & band)
{
    return band._upperDiagonal;
}

// ----------------------------------------------------------------------------
// Function getVerticalPos()
// ----------------------------------------------------------------------------

template <typename TBandSpec>
inline typename Position<Band<TBandSpec> >::Type
getLowerDiagonal(Band<TBandSpec> const & band)
{
    return band._lowerDiagonal;
}

// ----------------------------------------------------------------------------
// Function setHorizontalPos()
// ----------------------------------------------------------------------------

template <typename TBandSpec>
inline void
setUpperDiagonal(Band<BandSwitchedOn<TBandSpec> > & band,
                 typename Position<Band<TBandSpec> >::Type const & upperDiagonal)
{
    band._upperDiagonal = upperDiagonal;
    SEQAN_ASSERT(_checkValidBand(band));
}

// ----------------------------------------------------------------------------
// Function setVerticalPos()
// ----------------------------------------------------------------------------

template <typename TBandSpec>
inline void
setLowerDiagonal(Band<BandSwitchedOn<TBandSpec> > & band,
               typename Position<Band<TBandSpec> >::Type const & lowerDiagonal)
{
    band._lowerDiagonal = lowerDiagonal;
    SEQAN_ASSERT(_checkValidBand(band));
}

// ----------------------------------------------------------------------------
// Function _checkValidBand()
// ----------------------------------------------------------------------------

// This function checks whether the upper diagonal of the band is right of the
// lower diagonal.
// The function returns true if the horizontal position (upper diagonal) is
// greater or equal than the vertical position (lower diagonal).

template <typename TBandSpec>
inline bool _checkValidBand(Band<TBandSpec> const & band)
{
    return band._upperDiagonal >= band._lowerDiagonal;
}

// ----------------------------------------------------------------------------
// Function _isBandEnabled()
// ----------------------------------------------------------------------------

template <typename TBandSpec>
inline bool
_isBandEnabled(Band<BandSwitchedOn<TBandSpec> > const & /*band*/)
{
    return true;
}

inline bool
_isBandEnabled(Band<BandSwitchedOff > const & /*band*/)
{
    return false;
}

// ----------------------------------------------------------------------------
// Function _isWideBand()
// ----------------------------------------------------------------------------

template <typename TBandSpec, typename TSequenceH, typename TSequenceV>
inline bool
_isWideBand(Band<BandSwitchedOn<TBandSpec> > const & band,
            TSequenceH const & /*seqH*/,
            TSequenceV const & seqV)
{
    typedef typename MakeSigned<typename Size<TSequenceH>::Type>::Type TSignedSize;
    if (static_cast<TSignedSize>(length(seqV) + getLowerDiagonal(band)) < getUpperDiagonal(band))
        return true;
    return false;
}

template <typename TSequenceH, typename TSequenceV>
inline bool
_isWideBand(Band<BandSwitchedOff> const & /*band*/,
            TSequenceH const & /*seqH*/,
            TSequenceV const & /*seqV*/)
{
    return false;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_BAND_H_
