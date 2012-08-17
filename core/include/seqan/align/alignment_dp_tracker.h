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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_TRACKER_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_TRACKER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct ScoreOnly_{};
typedef Tag<ScoreOnly_> ScoreOnly;

struct RowMax_{};
typedef Tag<RowMax_> RowMax;

struct GridCoordinate
{
    size_t     x;
    size_t     y;

    GridCoordinate() : x(0), y(0) {}

    template <typename TPosition>
    GridCoordinate(TPosition const & posDim0, TPosition const & posDim1)
        :   x(posDim0),
            y(posDim1)
    {}
};

template<typename TScoreValue, typename TTraceIter, typename TSpec>
struct Tracker
{};

template<typename TScoreValue, typename TTraceIter>
struct Tracker<TScoreValue, TTraceIter, ScoreOnly >
{
    typedef String<TTraceIter> TTracePoints;

    TScoreValue   _maxScore;
    TTracePoints  _tracePoints;

    Tracker() : _maxScore(+MinValue<TScoreValue>::VALUE), _tracePoints()
    {
    }

    inline void operator()(TScoreValue const & score, TTraceIter const & /*tracePoint*/)
    {
        if (score > _maxScore)
        {
                // std::cout << "Debugging: " << _maxScore << " to " << dpValue << std::endl;
            _maxScore = score;
        }
    }
};

template<typename TScoreValue, typename TTraceIter>
struct Tracker<TScoreValue, TTraceIter, Default>
{
    typedef String<TTraceIter> TTracePoints;

    TScoreValue _maxScore;
    TTracePoints _tracePoints;

    Tracker() : _maxScore(+MinValue<TScoreValue>::VALUE), _tracePoints()
    {
        resize(_tracePoints, 1, Exact());
    }

    inline void operator()(TScoreValue const & score, TTraceIter const & tracePoint)
    {
            //new maximum found reset the trace point
        if (score > _maxScore)
        {
            _maxScore = score;
            _tracePoints[0] = tracePoint;
        }
    }
};

template<typename TScoreValue, typename TTraceIter>
    struct Tracker<TScoreValue, TTraceIter, RowMax > : public GridCoordinate
{
    size_t _minColNum;
    String<Pair<TScoreValue, size_t> > _maxScorePerRow;

    Tracker() : _minColNum(0), _maxScorePerRow()
    {
    }


    template <typename TRowSize>
    Tracker(size_t const & minColNum, TRowSize const & rowSize) : _minColNum(minColNum), _maxScorePerRow()
    {
        resize(_maxScorePerRow, rowSize, 0, Exact());
    }

    template<typename TIndex>
    inline void operator()(TScoreValue const & score, TTraceIter const & /*traceIter*/)
    {
            //new maximum found reset the trace point
        // TODO (rmaerker): change the way to compute the max row
        if (_minColNum <= x)
        {
            if (score > _maxScorePerRow[y].i1)
            {
                _maxScorePerRow[y].i1 = score;
                _maxScorePerRow[y].i2 = x;
            }
        }
        ++y;
    }
};

template<typename TScoreValue, typename TTraceIter>
struct Tracker<TScoreValue, TTraceIter, WatermanEggert>
{
    typedef String<TTraceIter> TTracePoints;

    TScoreValue _maxScore;
    TTracePoints _tracePoints;

    Tracker() :
    _maxScore(std::numeric_limits<TScoreValue>::min()), _tracePoints()
    {
    }

    inline void operator()(TScoreValue const & dpValue, TTraceIter const & tracePoint)
    {
            //new maximum found, replace old trace points with current
            //if same maximum found add new trace point
        if (dpValue == _maxScore)
        {
            appendValue(_tracePoints, tracePoint);
        }
        if (dpValue > _maxScore)
        {
            _maxScore = dpValue;
            clear(_tracePoints);
            appendValue(_tracePoints, tracePoint);
        }
    }
};


template<typename TScoreValue, typename TTraceIter>
struct Tracker<TScoreValue, TTraceIter, BandedChain > : public GridCoordinate
{
    typedef String<TTraceIter> TTracePoints;
    typedef DPValue<TScoreValue, AffineGaps> TDPValue;
    typedef String<TDPValue> TInitString;

    TScoreValue _maxScore;      // the maximal score detected
    size_t _posDim0;            // defines the column at which the values for _dim1Init should be tracked
    size_t _posDim1;            // defines the row at which the values for _dim0Init should be tracked
    TInitString _initDim0Curr;     // the horizontal initialization array of the current matrix
    TInitString _initDim1Curr;     // the vertical initialization array of the current matrix
    TInitString _initDim0Next;     // the horizontal initialization array of the next matrix
    TInitString _initDim1Next;     // the vertical initialization array of the next matrix
    TTracePoints _tracePoints;    // the array containing all tracepoints having the maximal score


    Tracker() : GridCoordinate(),
                _maxScore(MinValue<TScoreValue>::VALUE),
                _posDim0(0),
                _posDim1(0),
                _initDim0Curr(),
                _initDim1Curr(),
                _initDim0Next(),
                _initDim1Next(),
                _tracePoints()
    {
    }

    template <typename TSize, typename TPosition>
    Tracker(TSize const & lengthDim0Curr,
            TSize const & lengthDim1Curr,
            TPosition const & posDim0,
            TPosition const & posDim1,
            TSize const & lengthDim0Next,
            TSize const & lengthDim1Next) :
                GridCoordinate(),
                _maxScore(MinValue<TScoreValue>::VALUE),
                _posDim0(),
                _posDim1(),
                _initDim0Curr(),
                _initDim1Curr(),
                _initDim0Next(),
                _initDim1Next(),
                _tracePoints()
    {
        _init(*this, lengthDim0Curr, lengthDim1Curr, posDim0, posDim1, lengthDim0Next, lengthDim1Next);
    }

    inline void operator()(TScoreValue const & score, TTraceIter const & tracePoint)
    {
        // only track cells if they are within the correct region
        if (x >= _posDim0)
        {
            //new maximum found, replace old trace points with current
            //if same maximum found add new trace point
            if (score == _maxScore)
            {
                appendValue(_tracePoints, tracePoint);
            }
            if (score > _maxScore)
            {
                _maxScore = score;
                resize(_tracePoints,0);    // clear current trace points
                appendValue(_tracePoints, tracePoint); // add new value.
            }
        }
    }

    // this part fills the horizontal and vertical initialization vectors for next matrix in chain
    inline void operator()(TDPValue const & dpValue)
    {
        if (x == _posDim0)
        { // is in correct column
            if (y >= _posDim1)
            { // is within correct range in dim1
                _initDim1Next[y - _posDim1] = dpValue;
            }
        }
        else if (x >= _posDim0)
        {
            if (y == _posDim1)
            { // is in correct row
                _initDim0Next[x - _posDim0] = dpValue;
            }
        }
        ++y; // move further down in grid
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TAlignmentProfile>
struct GetTrackerSpec
{
    typedef Nothing Type;
};

template <typename TAlignCore, typename TGapSpec>
struct GetTrackerSpec<AlignmentProfile<TAlignCore,  TGapSpec, TracebackSwitchedOff> >
{
    typedef ScoreOnly Type;
};

template <typename TAlignCore, typename TGapSpec>
struct GetTrackerSpec<AlignmentProfile<TAlignCore, TGapSpec, TracebackSwitchedOn> >
{
    typedef Default Type;
};

template <typename TGapSpec>
struct GetTrackerSpec<AlignmentProfile<Local<WatermanEggert>,  TGapSpec, TracebackSwitchedOn> >
{
    typedef WatermanEggert Type;
};

template <typename TGapSpec>
struct GetTrackerSpec<AlignmentProfile<Global<RowMax>,  TGapSpec, TracebackSwitchedOn> >
{
    typedef RowMax Type;
};

template <typename TGapSpec>
struct GetTrackerSpec<AlignmentProfile<Global<BandedChain>,  TGapSpec, TracebackSwitchedOn> >
{
    typedef BandedChain Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// function _init()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TIter, typename TSize1, typename TPosition, typename TSize2>
inline void
_init(Tracker<TScoreValue, TIter, BandedChain > & tracker,
        TSize1 const & lengthDim0Curr,
        TSize1 const & lengthDim1Curr,
        TPosition const & posDim0,
        TPosition const & posDim1,
        TSize2 const & lengthDim0Next,
        TSize2 const & lengthDim1Next)
{
    typedef Tracker<TScoreValue, TIter, BandedChain> TTracker;
    typedef typename TTracker::TDPValue TDPValue;

    tracker.x = 0;
    tracker.y = 0;
    tracker._posDim0 = posDim0;
    tracker._posDim1 = posDim1;
    // initialize the current initialization values of the tracker
    arrayFill(begin(tracker._initDim0Curr), end(tracker._initDim0Curr), TDPValue(true));
    arrayFill(begin(tracker._initDim1Curr), end(tracker._initDim1Curr), TDPValue(true));
    arrayFill(begin(tracker._initDim0Next), end(tracker._initDim0Next), TDPValue(true));
    arrayFill(begin(tracker._initDim1Next), end(tracker._initDim1Next), TDPValue(true));

    // check if the value needs to be resized ... (can be longer but not smaller)
    if ((TSize1) length(tracker._initDim0Curr) < lengthDim0Curr)
    {
        resize(tracker._initDim0Curr, lengthDim0Curr, TDPValue(true), Exact());
    }
    if ((TSize1) length(tracker._initDim1Curr) < lengthDim1Curr)
    {
        resize(tracker._initDim1Curr, lengthDim1Curr, TDPValue(true), Exact());
    }
    if ((TSize2) length(tracker._initDim0Next) <= lengthDim0Next)
    {
        resize(tracker._initDim0Next, lengthDim0Next, TDPValue(true), Exact());
    }
    if ((TSize2) length(tracker._initDim1Next) <= lengthDim1Next)
    {
        resize(tracker._initDim1Next, lengthDim1Next, TDPValue(true), Exact());
    }
}

// functions for the GridCursor
// ----------------------------------------------------------------------------
// function goBeginNextColumn()
// ----------------------------------------------------------------------------

template <typename TTracker, typename TPosition, typename TBand>
inline void
goBeginNextColumn(TTracker const & /*cursor*/,
        TPosition const & /*size*/,
        TBand const & /*band*/)
{
    // nothing to do
}

template <typename TScoreValue, typename TIter, typename TPosition, typename TBandSpec>
inline void
goBeginNextColumn(Tracker<TScoreValue, TIter, BandedChain> & tracker,
        TPosition const & pos,
        Band<BandSwitchedOn<TBandSpec> > const &/*band*/)
{   // in banded version the span of begin position of vertical sequence is walking down the grid
    // the further the band goes right in horizontal direction
    ++tracker.x;
    tracker.y = pos + 1;
}

template <typename TScoreValue, typename TIter, typename TPosition>
inline void
goBeginNextColumn(Tracker<TScoreValue, TIter, BandedChain> & tracker,
        TPosition const & /*pos*/,
        Band<BandSwitchedOff> const & /*band*/ )
{
    ++tracker.x;
    tracker.y = 0;
}

template <typename TScoreValue, typename TIter, typename TPosition>
inline void
goBeginNextColumn(Tracker<TScoreValue, TIter, RowMax> & tracker, TPosition const & pos)
{
    ++tracker.x;
    tracker.y = pos + 1;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_TRACKER_H_
