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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_VALUE_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_VALUE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ForbiddenDirection
// ----------------------------------------------------------------------------

struct ForbiddenDirection
{
    static const uint8_t NONE = 0;  // no direction is forbidden
    static const uint8_t HORIZONTAL = 1; // horizontal cell is forbidden
    static const uint8_t VERTICAL = 2; // vertical cell is forbidden
};


template <typename TScoreValue, typename TConfig>
struct DPValue<TScoreValue, LinearGaps, TConfig> : TConfig
{
    TScoreValue _score;

    // default ctor
    DPValue() : TConfig(), _score(0) {}

    // copy ctor
    DPValue(DPValue<TScoreValue, LinearGaps, DPValueConfigDefault> const & other) :
            TConfig(), _score(other._score)
    {
    }

    // copy ctor for different configuration
    DPValue(DPValue<TScoreValue, LinearGaps, DPValueConfigForbidden> const & other) :
            TConfig(other._forbidden),
            _score(other._score)
    {
    }

    // ctor with score
    DPValue(TScoreValue const & score) :
            TConfig(),
            _score(score)
    {
    }

    // ctor with forbidden flag
    DPValue(bool const & forbidden) :
            TConfig(forbidden),
            _score(0)
    {
    }

    DPValue<TScoreValue, LinearGaps, DPValueConfigDefault> &
    operator=(DPValue<TScoreValue, LinearGaps, DPValueConfigDefault> const & other)
    {
        if (this != &other)
            _score = other._score;
        return *this;
    }

    DPValue<TScoreValue, LinearGaps, DPValueConfigForbidden> &
    operator=(DPValue<TScoreValue, LinearGaps, DPValueConfigForbidden> const & other)
    {
        if (this != &other)
        {
            this->_forbidden = other._forbidden;
            _score = other._score;
        }
        return *this;
    }

    DPValue &
    operator=(TScoreValue const & score)
    {
        _score = score;
        return *this;
    }
};

template <typename TScoreValue, typename TConfig>
struct DPValue<TScoreValue, AffineGaps, TConfig> : TConfig
{
    TScoreValue _score;
    TScoreValue _scoreHorizontal;
    TScoreValue _scoreVertical;
    uint8_t     _dirForbidden;


    DPValue() :  TConfig(), _score(0), _scoreHorizontal(0), _scoreVertical(0), _dirForbidden(ForbiddenDirection::NONE) {}


    DPValue(DPValue<TScoreValue, AffineGaps, DPValueConfigDefault> const & other) :
            TConfig(),
            _score(other._score),
            _scoreHorizontal(other._scoreHorizontal),
            _scoreVertical(other._scoreVertical),
            _dirForbidden(other._dirForbidden)
    {
    }

    DPValue(DPValue<TScoreValue, AffineGaps, DPValueConfigForbidden> const & other) :
            TConfig(other._forbidden),
            _score(other._score),
            _scoreHorizontal(other._scoreHorizontal),
            _scoreVertical(other._scoreVertical),
            _dirForbidden(other._dirForbidden)
    {
    }

    DPValue(TScoreValue const & score) :
            TConfig(),
            _score(score),
            _scoreHorizontal(0),
            _scoreVertical(0),
            _dirForbidden(ForbiddenDirection::NONE)
    {
    }

    DPValue(bool const & forbidden) :
            TConfig(forbidden),
            _score(0),
            _scoreHorizontal(0),
            _scoreVertical(0),
            _dirForbidden(ForbiddenDirection::NONE)
    {
    }

    DPValue<TScoreValue, AffineGaps, DPValueConfigDefault> &
    operator=(DPValue<TScoreValue, AffineGaps, DPValueConfigDefault> const & other)
    {
        if (*this != other)
        {
            _score = other._score;
            _scoreHorizontal = other._scoreHorizontal;
            _scoreVertical = other._scoreVertical;
            _dirForbidden = other._dirForbidden;
        }
        return *this;
    }

    DPValue<TScoreValue, AffineGaps, DPValueConfigForbidden> &
    operator=(DPValue<TScoreValue, AffineGaps, DPValueConfigForbidden> const & other)
    {
        if (*this != other)
        {
            this->_forbidden = other._forbidden;
            _score = other._score;
            _scoreHorizontal = other._scoreHorizontal;
            _scoreVertical = other._scoreVertical;
            _dirForbidden = other._dirForbidden;
        }
        return *this;
    }

    DPValue &
    operator=(TScoreValue const & score)
    {
        _score = score;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Value
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TDPValueSpec, typename TDPValueConfig>
struct Value<DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> >
{
    typedef TScoreValue Type;
};

template <typename TScoreValue, typename TDPValueSpec, typename TDPValueConfig>
struct Value<DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> const >
{
    typedef TScoreValue const Type;
};

// ----------------------------------------------------------------------------
// Reference
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TDPValueSpec, typename TDPValueConfig>
struct Reference<DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> >
{
    typedef TScoreValue & Type;
};

template <typename TScoreValue, typename TDPValueSpec, typename TDPValueConfig>
struct Reference<DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> const >
{
    typedef TScoreValue const & Type;
};

// ----------------------------------------------------------------------------
// Spec
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TDPValueSpec, typename TDPValueConfig>
struct Spec<DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> >
{
    typedef TDPValueSpec Type;
};

template <typename TScoreValue, typename TDPValueSpec, typename TDPValueConfig>
struct Spec<DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> const> :
            Spec<DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> > {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function operator<<()
// ----------------------------------------------------------------------------

template <typename TStream, typename TScoreValue, typename TDPValueSpec, typename TDPValueConfig>
inline TStream &
operator<<(TStream & stream, DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> const & value)
{
    if (isForbidden(value))
    {
        stream << "X";
    }
    else
    {
        stream << getScore(value);
    }
    return stream;
}

template <typename TStream, typename TScoreValue, typename TDPValueConfig>
inline TStream &
operator<<(TStream & stream, DPValue<TScoreValue, AffineGaps, TDPValueConfig> const & value)
{
    if (isForbidden(value))
    {
        stream << "X";
    }
    else
    {
        stream << getScore(value) << "(" << value._scoreHorizontal << ","<< value._scoreVertical <<")";
    }
    return stream;
}


// ----------------------------------------------------------------------------
// function operator==()                                           [LinearGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline bool
operator==(DPValue<TScoreValue, LinearGaps, DPValueConfigDefault> const & left,
           DPValue<TScoreValue, LinearGaps, DPValueConfigDefault> const & right)
{
    return left._score == right._score;
}

template <typename TScoreValue>
inline bool
operator==(DPValue<TScoreValue, LinearGaps, DPValueConfigForbidden> const & left,
           DPValue<TScoreValue, LinearGaps, DPValueConfigForbidden> const & right)
{
    if (left._forbidden != right._forbidden)
    {
        return false;
    }
    return left._score == right._score;
}

// ----------------------------------------------------------------------------
// function operator==()                                           [AffineGaps]
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline bool
operator==(DPValue<TScoreValue, AffineGaps, DPValueConfigForbidden> const & left,
           DPValue<TScoreValue, AffineGaps, DPValueConfigForbidden> const & right)
{
    if (left._forbidden != right._forbidden)
    {
        return false;
    }
    if (left._score != right._score)
    {
        return false;
    }
    if (left._scoreHorizontal != right._scoreHorizontal)
    {
        return false;
    }
    if (left._scoreVertical != right._scoreVertical)
    {
        return false;
    }
    if (left._dirForbidden != right._dirForbidden)
    {
        return false;
    }
    return true;
}

template <typename TScoreValue>
inline bool
operator==(DPValue<TScoreValue, AffineGaps, DPValueConfigDefault> const & left,
           DPValue<TScoreValue, AffineGaps, DPValueConfigDefault> const & right)
{
    if (left._score != right._score)
    {
        return false;
    }
    if (left._scoreHorizontal != right._scoreHorizontal)
    {
        return false;
    }
    if (left._scoreVertical != right._scoreVertical)
    {
        return false;
    }
    if (left._dirForbidden != right._dirForbidden)
    {
        return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// function operator!=()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGaps, typename TConfig>
inline bool
operator!=(DPValue<TScoreValue, TGaps, TConfig> const & left,
           DPValue<TScoreValue, TGaps, TConfig> const & right)
{
    return !(left == right);
}

// ----------------------------------------------------------------------------
// function getScore()
// ----------------------------------------------------------------------------

template<typename TScoreValue, typename TDPValueSpec, typename TDPValueConfig>
inline TScoreValue const &
getScore(DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> const & value)
{
    return value._score;
}

// ----------------------------------------------------------------------------
// function setScore()
// ----------------------------------------------------------------------------

template<typename TScoreValue, typename TDPValueSpec, typename TDPValueConfig>
inline TScoreValue const &
setScore(DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> & value, TScoreValue const & score)
{
    return value._score = score;
}

// ----------------------------------------------------------------------------
// function isForbidden()
// ----------------------------------------------------------------------------

template<typename TScoreValue, typename TDPValueSpec>
inline bool
isForbidden(DPValue<TScoreValue, TDPValueSpec, DPValueConfigDefault> const & /*value*/)
{
    return false;
}

template<typename TScoreValue, typename TDPValueSpec>
inline bool
isForbidden(DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> const & value)
{
    return value._forbidden;
}

// ----------------------------------------------------------------------------
// function setForbidden()
// ----------------------------------------------------------------------------

template<typename TScoreValue, typename TDPValueSpec>
inline void
setForbidden(DPValue<TScoreValue, TDPValueSpec, DPValueConfigDefault> & /*value*/,
             typename TraceBitMask::Type const & /*origin*/)
{
    // nothing to do
}

template<typename TScoreValue, typename TDPValueSpec>
inline void
setForbidden(DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> & value,
             bool const & forbidden)
{
    value._forbidden = forbidden;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_VALUE_H_
