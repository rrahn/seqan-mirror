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

// ============================================================================
// Functions
// ============================================================================

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

template <typename TScoreValue, typename TIter, typename TPosition>
inline void
goBeginNextColumn(Tracker<TScoreValue, TIter, RowMax> & tracker, TPosition const & pos)
{
    ++tracker.x;
    tracker.y = pos + 1;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_TRACKER_H_
