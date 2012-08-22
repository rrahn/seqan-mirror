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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_TRACESEGMENT_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_TRACESEGMENT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TPosition, typename TSize>
class TraceSegment
{
public:
    typedef typename TraceBitMask::Type TTraceValue;

    TPosition _horizontalBeginPos;
    TPosition _verticalBeginPos;
    TSize _length;
    TTraceValue _traceValue;

    TraceSegment() : _horizontalBeginPos(0), _verticalBeginPos(0), _length(0), _traceValue(+TraceBitMask::NONE){}

    TraceSegment(TraceSegment const & other) : _horizontalBeginPos(other._horizontalBeginPos),
    _verticalBeginPos(other._verticalBeginPos),
    _length(other._length),
    _traceValue(other._traceValue) {}

    TraceSegment(TPosition const & horizontalBeginPos, TPosition const & verticalBeginPos, TSize const & length,
                 TTraceValue const & traceValue) : _horizontalBeginPos(horizontalBeginPos),
    _verticalBeginPos(verticalBeginPos),
    _length(length),
    _traceValue(traceValue) {}

    TraceSegment &
    operator=(TraceSegment const & other)
    {
        if (*this != other)
        {
            _horizontalBeginPos = other._horizontalBeginPos;
            _verticalBeginPos = other._verticalBeginPos;
            _length = other._length;
            _traceValue = other._traceValue;
        }
        return *this;
    }



};


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Position
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
struct Position<TraceSegment<TPosition, TSize> >
{
    typedef TPosition Type;
};

template <typename TPosition, typename TSize>
struct Position<TraceSegment<TPosition, TSize> const> :
    Position<TraceSegment<TPosition, TSize> > {};

// ----------------------------------------------------------------------------
// Size
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
struct Size<TraceSegment<TPosition, TSize> >
{
    typedef TSize Type;
};

template <typename TPosition, typename TSize>
struct Size<TraceSegment<TPosition, TSize> const> :
    Size<TraceSegment<TPosition, TSize> > {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// function getBeginHorizontal()
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
inline TPosition
getBeginHorizontal(TraceSegment<TPosition, TSize> const & traceSegment)
{
    return traceSegment._horizontalBeginPos;
}

// ----------------------------------------------------------------------------
// function getBeginVertical()
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
inline TPosition
getBeginVertical(TraceSegment<TPosition, TSize> const & traceSegment)
{
    return traceSegment._verticalBeginPos;
}

// ----------------------------------------------------------------------------
// function getEndHorizontal()
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
inline TPosition
getEndHorizontal(TraceSegment<TPosition, TSize> const & traceSegment)
{
    if (traceSegment._traceValue & (TraceBitMask::HORIZONTAL | TraceBitMask::DIAGONAL ))
    {
        return traceSegment._horizontalBeginPos + traceSegment._length;
    }
    return traceSegment._horizontalBeginPos;
}

// ----------------------------------------------------------------------------
// function getEndVertical()
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
inline TPosition
getEndVertical(TraceSegment<TPosition, TSize> const & traceSegment)
{
    if (traceSegment._traceValue & (TraceBitMask::VERTICAL | TraceBitMask::DIAGONAL))
    {
        return traceSegment._verticalBeginPos + traceSegment._length;
    }
    return traceSegment._verticalBeginPos;
}

// ----------------------------------------------------------------------------
// function length()
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
inline TSize
length(TraceSegment<TPosition, TSize> const & traceSegment)
{
    return traceSegment._length;
}

// TODO (rmaerker): remove debug code

template <typename TTraceValue>
String<char> translateTraceValue(TTraceValue const & traceValue)
{
    String<char> transcript;

    if ((traceValue & TraceBitMask::DIAGONAL) == TraceBitMask::DIAGONAL)
    {
        append(transcript, 'D');
    }
    if ((traceValue & TraceBitMask::VERTICAL) == TraceBitMask::VERTICAL)
    {
        append(transcript, 'V');
    }
    if ((traceValue & TraceBitMask::HORIZONTAL) == TraceBitMask::HORIZONTAL)
    {
        append(transcript, 'H');
    }
    if ((traceValue & TraceBitMask::VERTICAL_OPEN) == TraceBitMask::VERTICAL_OPEN)
    {
        append(transcript, 'v');
    }
    if ((traceValue & TraceBitMask::HORIZONTAL_OPEN) == TraceBitMask::HORIZONTAL_OPEN)
    {
        append(transcript, 'h');
    }
    if ((traceValue) == TraceBitMask::NONE)
    {
        append(transcript, '0');
    }
    return transcript;
}

template <typename TStream, typename TSize, typename TPosition>
TStream & operator<<(TStream & stream, TraceSegment<TSize, TPosition> const & traceSegment)
{
    stream << translateTraceValue(traceSegment._traceValue) << "-";
    stream << "(" << traceSegment._horizontalBeginPos << ", " << traceSegment._verticalBeginPos << ", " <<
                traceSegment._length << ")";
    return stream;
}

template <typename TPosition, typename TSize>
inline bool operator==(TraceSegment<TPosition, TSize> const & left, TraceSegment<TPosition, TSize> const & right)
{
    if (left._horizontalBeginPos != right._horizontalBeginPos)
    {
        return false;
    }
    if (left._verticalBeginPos != right._verticalBeginPos)
    {
        return false;
    }
    if (left._length != right._length)
    {
        return false;
    }
    if (left._traceValue != right._traceValue)
    {
        return false;
    }
    return true;
}

template <typename TPosition, typename TSize>
inline bool operator!=(TraceSegment<TPosition, TSize> const & left, TraceSegment<TPosition, TSize> const & right)
{
    return !(left == right);
}

//---------------------------------------------------
// function recordSegment()
//---------------------------------------------------

template <typename TTraceSegments, typename TPosition, typename TSize, typename TTraceValue>
inline void recordSegment(TTraceSegments & traceSegments,
                          TPosition const & horizontalBeginPos,
                          TPosition const & verticalBeginPos,
                          TSize const & segmentLength,
                          TTraceValue const & traceValue)
{
    typedef typename Value<TTraceSegments>::Type TTraceSegment;
    if (segmentLength == 0)
    {
        return; // we don't store empty segments
    }
    if (traceValue & TraceBitMask::DIAGONAL)
    {
        appendValue(traceSegments, TTraceSegment(horizontalBeginPos, verticalBeginPos, segmentLength, +TraceBitMask::DIAGONAL));
    }
    else if (traceValue & TraceBitMask::VERTICAL)
    {
        appendValue(traceSegments, TTraceSegment(horizontalBeginPos, verticalBeginPos, segmentLength, +TraceBitMask::VERTICAL));
    }
    else if (traceValue & TraceBitMask::HORIZONTAL)
    {
        appendValue(traceSegments, TTraceSegment(horizontalBeginPos, verticalBeginPos, segmentLength, +TraceBitMask::HORIZONTAL));
    }
    // everything else is not tracked.
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_TRACESEGMENT_H_
