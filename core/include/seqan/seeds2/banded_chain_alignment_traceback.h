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

#ifndef CORE_INCLUDE_SEQAN_SEEDS2_BANDED_CHAIN_ALIGNMENT_TRACEBACK_H_
#define CORE_INCLUDE_SEQAN_SEEDS2_BANDED_CHAIN_ALIGNMENT_TRACEBACK_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TInitSpotArray>
inline bool
_existInitSpot(TInitSpotArray const & initSpots,
        typename Value<TInitSpotArray>::Type const & elem)
{
    typedef typename Iterator<TInitSpotArray const>::Type TInitSpotArrayIterator;

    for (TInitSpotArrayIterator it = begin(initSpots); it != end(initSpots); ++it)
    {
        if (value(it) == elem)
        {
            return true;
        }
    }
    return false;
}

template <typename TTraceSet, typename TGridPoint>
inline void
_adaptLocalTracesToGlobalGrid(TTraceSet & traceSet, TGridPoint const & gridBegin)
{
    typedef typename Value<TTraceSet>::Type TTraceString;
    typedef typename Iterator<TTraceString>::Type TTraceStringIterator;

    for (unsigned i = 0; i < length(traceSet); ++i)
    {
        for (TTraceStringIterator it = begin(traceSet[i]); it != end(traceSet[i]); ++it)
        {
            value(it)._horizontalBeginPos += gridBegin.i1;
            value(it)._verticalBeginPos += gridBegin.i2;
        }
    }
}

template <typename TTraceSet>
inline void glueTraces(TTraceSet & globalTraces, TTraceSet const & localTraces)
{
    typedef typename Value<TTraceSet>::Type TTraceSegments;
    typedef typename Value<TTraceSegments>::Type TTraceSegment;
    typedef typename Size<TTraceSegments>::Type TSize;

    typedef typename Iterator<TTraceSegments const>::Type TTraceSegmentsConstIterator;
    typedef typename Iterator<TTraceSegments>::Type TTraceSegmentsIterator;

    bool isGlued = false;

    size_t lengthGlobalTraces = length(globalTraces);
    String<unsigned> elementsToErase;

    for (unsigned j = 0; j < lengthGlobalTraces; ++j)
    {
        // traceback from back to front -> first trace segment is at end of sequences
        TTraceSegment globalTraceEndPoint = value(begin(globalTraces[j]));
        TSize numOfCurrElements = length(globalTraces[j]);
        // check for all existing paths if traces can be glued together.
        bool isConnected = false;
        for (unsigned i = 0; i < length(localTraces); ++i)
        {
            // traceback from back to front -> last trace segment is at beginning of sequences.
            TTraceSegment localTraceBeginPoint = value(end(localTraces[i]) - 1);
            // trace segments match in horizontal position
            if (getEndHorizontal(globalTraceEndPoint) == getBeginHorizontal(localTraceBeginPoint))
            {
                // trace segments match in vertical position
                if (getEndVertical(globalTraceEndPoint) == getBeginVertical(localTraceBeginPoint))
                { // found a glue point between local and global trace
                    TSize numOfElementsToAdd = length(localTraces[i]);
                    if (isConnected)
                    {
                        // create a new traceback. track
                        TTraceSegments newTraceTrack;
                        resize(newTraceTrack, numOfCurrElements + numOfElementsToAdd, Generous());
                        arrayMoveForward(begin(localTraces[i]), end(localTraces[i]), begin(newTraceTrack));
                        arrayCopyForward(end(globalTraces[j]) - numOfCurrElements, end(globalTraces[j]), begin(newTraceTrack) + numOfElementsToAdd);
                        appendValue(globalTraces, newTraceTrack);
                        continue;
                    }
                    // resize global traces such that new elements fit into it
                    resize(globalTraces[j], numOfCurrElements + numOfElementsToAdd, Generous());

                    // shift old values to correct position in array after new elements will be added.
                    // use backward move in order to avoid overwriting when ranges overlap
                    arrayMoveBackward(begin(globalTraces[j]), begin(globalTraces[j]) + numOfCurrElements,
                                      begin(globalTraces[j]) + numOfElementsToAdd);
                    // actually move new elements to global traces at it's begin
                    arrayMoveForward(begin(localTraces[i]), end(localTraces[i]), begin(globalTraces[j]));
                    isGlued = true;
                    isConnected = true;
                }
            }
        }
        if (!isConnected)
        {
            appendValue(elementsToErase, j);
        }
    }

    for (unsigned i = length(elementsToErase); i > 0; --i)
    {
        erase(globalTraces, elementsToErase[i-1]); // erase from behind to avoid accessing an element beyond the scope
    }
    SEQAN_ASSERT_EQ_MSG(isGlued, true, "Fatal error while trying to connect trace backs: No glue point available!");
}

template <typename TTarget, typename TInitSpot, typename TTraceIter, typename TTracker, typename TTraceMatrix,
        typename TSeq0, typename TSeq1,  typename TBand>
inline void
followTrace(TTarget & target,
            TInitSpot & initSpot,
            TTraceIter & traceIter,
            TTracker const & tracker,
            TTraceMatrix const & traceMatrix,
            TSeq0 const & seq0,
            TSeq1 const & seq1,
            TBand const & band,
            AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn> const & /*profile*/)
{
    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn> TProfile;
    typedef typename GetGapSpec<TProfile>::Type TGapSpec;
    typedef typename Position<TSeq0>::Type TSeq0Pos;
    typedef typename Position<TSeq1>::Type TSeq1Pos;
    typedef typename Infix<TSeq0>::Type TInfix0;
    typedef typename Infix<TSeq1>::Type TInfix1;
    typedef typename Difference<TTraceIter>::Type TDiff;

    // first follow the traceback up to the row and or column that is used for intialization for the consecutive
    // matrix.
    TDiff tracePositon = traceIter - static_cast<TTraceIter>(begin(traceMatrix));
    TSeq0Pos suffix0Pos = _determineSeq0Pos(tracePositon, seq0, seq1, band) - tracker._posDim0;
    TSeq1Pos suffix1Pos = _determineSeq1Pos(tracePositon, seq0, seq1, band) - tracker._posDim1;

    // TODO(rmaerker): we need the first and the last anchor here. but why...

    // Don't need the actual path as for this part it is determined during the traceback of the next matrix
    TTarget tmpTarget;
    _doFollowTrace(tmpTarget, traceIter, suffix0Pos, suffix1Pos, getColumnSize(seq0, seq1, band), band, TGapSpec());

    // store the initSpot that was detected in first part
    tracePositon = traceIter - static_cast<TTraceIter>(begin(traceMatrix));
    TSeq0Pos seq0Pos = _determineSeq0Pos(tracePositon, seq0, seq1, band);
    TSeq1Pos seq1Pos = _determineSeq1Pos(tracePositon, seq0, seq1, band);
    initSpot.i1 = seq0Pos - tracker._posDim0;
    initSpot.i2 = seq1Pos - tracker._posDim1;
    if (initSpot.i1 == 0)
    {
        initSpot.i3 = tracker._initDim1Next[initSpot.i2];
    }
    else if (initSpot.i2 == 0)
    {
        initSpot.i3 = tracker._initDim0Next[initSpot.i1];
    }
    // continue traceback to the end and actually store it
    _doFollowTrace(target, traceIter, seq0Pos, seq1Pos, getColumnSize(seq0, seq1, band), band, TGapSpec());
}

template <typename TTarget, typename TTraceIter, typename TTraceMatrix,
        typename TSeq0, typename TSeq1,  typename TBand>
inline void
followTrace(TTarget & target,
            TTraceIter & traceIter,
            TTraceMatrix const & traceMatrix,
            TSeq0 const & seq0,
            TSeq1 const & seq1,
            TBand const & band,
            AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn> const & /*profile*/)
{
    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn> TProfile;
    typedef typename GetGapSpec<TProfile>::Type TGapSpec;
    typedef typename Position<TSeq0>::Type TSeq0Pos;
    typedef typename Position<TSeq1>::Type TSeq1Pos;
    typedef typename Infix<TSeq0>::Type TInfix0;
    typedef typename Infix<TSeq1>::Type TInfix1;
    typedef typename Difference<TTraceIter>::Type TDiff;

    // first follow the traceback up to the row and or column that is used for intialization for the consecutive
    // matrix.

    TDiff tracePosition = traceIter - static_cast<TTraceIter>(begin(traceMatrix));
    TSeq0Pos seq0Pos = _determineSeq0Pos(tracePosition, seq0, seq1, band);
    TSeq1Pos seq1Pos = _determineSeq1Pos(tracePosition, seq0, seq1, band);

    // record overlaps at tail of sequences if any
    if (seq0Pos != length(seq0))
    {
        recordSegment(target, seq0Pos, seq1Pos, length(seq0) - seq0Pos, +TraceBitMask::HORIZONTAL);
    }
    if (seq1Pos != length(seq1))
    {
        recordSegment(target, seq0Pos, seq1Pos, length(seq1) - seq1Pos, +TraceBitMask::VERTICAL);
    }

    // actually follow the trace
    _doFollowTrace(target, traceIter, seq0Pos, seq1Pos, getColumnSize(seq0, seq1, band), band, TGapSpec());
}

template <typename TTraceSegSet, typename TInitSpots, typename TTracker, typename TTraceMatrix, typename TSeq0, typename TSeq1,
          typename TBand, typename TGapSpec, typename TTraceSpec>
inline void
computeTraceback(TTraceSegSet & localTraceSet,
                 TInitSpots & initSpots,
                 TTracker & tracker,
                 TTraceMatrix const & traceMatrix,
                 TSeq0 const & seq0,
                 TSeq1 const & seq1,
                 TBand const & band,
                 AlignmentProfile<Global<BandedChain>, TGapSpec, TTraceSpec> const & profile)
{
    typedef typename Infix<TSeq0>::Type TInfix0;
    typedef typename Infix<TSeq1>::Type TInfix1;
    typedef typename TTracker::TTracePoints TTracePoints;
    typedef typename Iterator<TTracePoints>::Type TTracePointIterator;
    typedef typename Value<TInitSpots>::Type TInitSpot;
    typedef typename Value<TTraceSegSet>::Type TTraceSegments;

    clear(initSpots);
    clear(localTraceSet);

    TTracePointIterator it = begin(tracker._tracePoints);
    TTracePointIterator itEnd = end(tracker._tracePoints);

    for (; it != itEnd; ++it)
    {
        TTraceSegments tmpTraceSegments;
        TInitSpot tmpInitSpot;
        followTrace(tmpTraceSegments, tmpInitSpot, value(it), tracker, traceMatrix, seq0, seq1, band, profile);

        if (!_existInitSpot(initSpots, tmpInitSpot))
        {
            appendValue(localTraceSet, tmpTraceSegments, Exact());
            appendValue(initSpots, tmpInitSpot);
        }
    }
    SEQAN_ASSERT_EQ(length(initSpots), length(localTraceSet));
}

template <typename TTraceSegSet, typename TTracker, typename TTraceMatrix, typename TSeq0, typename TSeq1,
          typename TBand, typename TGapSpec, typename TTraceSpec>
inline void
computeTraceback(TTraceSegSet & localTraceSet,
                 TTracker & tracker,
                 TTraceMatrix const & traceMatrix,
                 TSeq0 const & seq0,
                 TSeq1 const & seq1,
                 TBand const & band,
                 AlignmentProfile<Global<BandedChain>, TGapSpec, TTraceSpec> const & profile)
{
    typedef typename Infix<TSeq0>::Type TInfix0;
    typedef typename Infix<TSeq1>::Type TInfix1;
    typedef typename TTracker::TTracePoints TTracePoints;
    typedef typename Iterator<TTracePoints>::Type TTracePointIterator;

    clear(localTraceSet);

    TTracePointIterator it = begin(tracker._tracePoints);
    TTracePointIterator itEnd = end(tracker._tracePoints);
    resize(localTraceSet, length(tracker._tracePoints), Exact());

    unsigned pos = 0;
    for (; it != itEnd; ++it, ++pos)
    {
        followTrace(localTraceSet[pos], value(it), traceMatrix, seq0, seq1, band, profile);
    }
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_SEEDS2_BANDED_CHAIN_ALIGNMENT_TRACEBACK_H_
