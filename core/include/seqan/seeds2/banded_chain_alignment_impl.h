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

#ifndef CORE_INCLUDE_SEQAN_SEEDS2_BANDED_CHAIN_ALIGNMENT_IMPL_H_
#define CORE_INCLUDE_SEQAN_SEEDS2_BANDED_CHAIN_ALIGNMENT_IMPL_H_

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

// Helper functions for banded chain alignment

// checks whether a given seeds crosses the first row defined in the grid spanned by the both given sequences.
template <typename TPosition, typename TSize>
inline bool
isCrossingBeginDim0(TPosition const & posDim0, TSize const & bandSize)
{
    typedef typename MakeSigned<TSize>::Type TSignedSize;

    if ((TSignedSize) posDim0 - (TSignedSize) bandSize <= 0)
    {
        return true;
    }
    return false;
}

// checks whether a given seeds crosses the first column defined in the grid spanned by the both given sequences.
template <typename TPosition, typename TSize>
inline bool
isCrossingBeginDim1(TPosition const & posDim1, TSize const & bandSize)
{
    typedef typename MakeSigned<TSize>::Type TSignedSize;

    if ((TSignedSize) posDim1 - (TSignedSize) bandSize <= 0)
    {
        return true;
    }
    return false;
}

// checks whether a given seeds crosses the last row defined in the grid spanned by the both given sequences.
template <typename TPosition, typename TSize, typename TSequence0>
inline bool
isCrossingEndDim0(TPosition const & posDim0, TSize const & bandSize, TSequence0 const & seq0)
{
    typedef typename MakeSigned<TSize>::Type TSignedSize;

    if ((TSignedSize) posDim0 + (TSignedSize) bandSize >= (TSignedSize) length(seq0))
    {
        return true;
    }
    return false;
}

// checks whether a given seeds crosses the last column defined in the grid spanned by the both given sequences.
template <typename TPosition, typename TSize, typename TSequence1>
inline bool
isCrossingEndDim1(TPosition const & posDim1, TSize const & bandSize, TSequence1 const & seq1)
{
    typedef typename MakeSigned<TSize>::Type TSignedSize;

    if ((TSignedSize) posDim1 + (TSignedSize) bandSize >= (TSignedSize) length(seq1))
    {
        return true;
    }
    return false;
}

// TODO (rmaerker): switch orientation of dim0/dim1 - from vertical/horizontal to horizontal/vertical in seeds module
template <typename TSeed>
inline typename Size<TSeed>::Type
getBandShiftHorizontal(TSeed const & seed)
{
    typedef typename Size<TSeed>::Type TSize;
    TSize bandSize = (getUpperDiagonal(seed) - (getBeginDim1(seed) - getBeginDim0(seed))); // must be greater or equal than zero
    SEQAN_ASSERT_GEQ(bandSize, TSize(0));
    return bandSize;
}

// TODO (rmaerker): switch orientation of dim0/dim1 - from vertical/horizontal to horizontal/vertical in seeds module
template <typename TSeed>
inline typename Size<TSeed>::Type
getBandShiftVertical(TSeed const & seed)
{
    typedef typename Size<TSeed>::Type TSize;
    TSize bandSize ((getBeginDim1(seed) - getBeginDim0(seed)) - getLowerDiagonal(seed)); // must be greater or equal than zero
    SEQAN_ASSERT_GEQ(bandSize, TSize(0));
    return bandSize;
}

// sets ths initialization spots discovered in the previous run over the predecessor matrix
template <typename TTracker, typename TInitSpots>
inline void
_setInitSpots(TTracker & tracker, TInitSpots const & initSpots)
{
    typedef typename Value<TInitSpots>::Type TInitSpot;

    for (unsigned i = 0; i < length(initSpots); ++i)
    {
        TInitSpot spot = initSpots[i];

        if (spot.i1 == 0)
        {
            tracker._initDim1Curr[spot.i2] = spot.i3;
        }
        else if (spot.i2 == 0)
        {
            tracker._initDim0Curr[spot.i1] = spot.i3;
        }
    }
}


template <typename TDPIter, typename TTracker, typename TBand, typename TColumnType>
inline void
_setInitValue(TDPIter & /*activeColIter*/,
               TTracker const & /*tracker*/,
               TBand const & /*band*/,
               TColumnType const /*columnType*/)
{
    // nothing to do
}

template <typename TDPIter, typename TScoreValue, typename TTraceIter>
inline void
_setInitValue(TDPIter & activeColIter,
               Tracker<TScoreValue, TTraceIter, BandedChain> const & tracker,
               Band<BandSwitchedOn<> > const & /*band*/,
               BandOpenColumn const & /*columnType*/)  // TODO(rmaerker): change to new values.
{
    value(activeColIter) = tracker._initDim0Curr[tracker.x];
}

template <typename TDPIter, typename TScoreValue, typename TTraceIter>
inline void
_setInitValue(TDPIter & activeColIter,
               Tracker<TScoreValue, TTraceIter, BandedChain> const & tracker,
               Band<BandSwitchedOff> const & /*band*/,
               FullColumn const /*columnType*/)
{
    value(activeColIter) = tracker._initDim0Curr[tracker.x];
}

template <typename TSpec, typename TConfig>
inline typename Position<Seed<TSpec, TConfig> >::Type
_getBeginDim0(Seed<TSpec, TConfig> const & seed)
{
    return getBeginDim0(seed) + 1;
}

template <typename TSpec, typename TConfig>
inline typename Position<Seed<TSpec, TConfig> >::Type
_getBeginDim1(Seed<TSpec, TConfig> const & seed)
{
    return getBeginDim1(seed) + 1;
}



template <typename TSeedSet, typename TBandSize>
inline typename Iterator<TSeedSet const, Standard>::Type
findFirstAnchor(TSeedSet const & seedSet, TBandSize const & bandWidth)
{
    // TODO (rmaerker): iterator should be annotated as BiDirectional Iterator see (std::set)
    typedef typename Iterator<TSeedSet const, Standard>::Type TIterator;
    typedef typename Value<TSeedSet const>::Type TSeed;

    SEQAN_ASSERT_GT_MSG(length(seedSet), 0u, "SeedSet is empty!");

    TIterator it = begin(seedSet);
    TIterator itEnd = end(seedSet);
    --itEnd;

    while(it != itEnd)
    {
        TSeed seed = value(++it);
        if (isCrossingBeginDim0(_getBeginDim0(seed), bandWidth))
        {
            continue;
        }
        else if (isCrossingBeginDim1(_getBeginDim1(seed), bandWidth))
        {
            continue;
        }
        else
        { // found seed which is not crossing the begin.
            SEQAN_ASSERT(it != static_cast<TIterator>(begin(seedSet)));
            return --it;
        }
    }
    return it;
}

template <typename TSeedSet, typename TSequence0, typename TSequence1, typename TBandSize>
inline typename Iterator<TSeedSet const, Standard>::Type
findLastAnchor(TSeedSet const & seedSet, TSequence0 const & seq0, TSequence1 const & seq1, TBandSize const & bandWidth)
{
    typedef typename Iterator<TSeedSet const, Standard>::Type TIterator;
    typedef typename Value<TSeedSet>::Type TSeed;

    SEQAN_ASSERT_GT_MSG(length(seedSet), 0u, "SeedSet is empty!");

    TIterator it = end(seedSet);
    --it;
    TIterator itEnd = begin(seedSet);

    while(it != itEnd)
    {
        TSeed seed = value(--it);
        if (isCrossingEndDim0(getEndDim0(seed), bandWidth, seq0))
        {
            continue;
        }
        else if (isCrossingEndDim1(getEndDim1(seed), bandWidth, seq1))
        {
            continue;
        }
        else
        { // found seed which is not crossing the end.
            return it;
        }
    }
    return it;
}

// Overloading functions to compute special alignment for banded chain algorithm


//---------------------------------------------------
// function track()
//---------------------------------------------------

// overloaded track function to implement different tracking method

template <typename TScoreValue, typename TTraceIter, typename TDPValue>
inline void track(Tracker<TScoreValue, TTraceIter, BandedChain> & tracker,
                  TDPValue const & dpValue,
                  TTraceIter const & traceIter,
                  True const & /*trackScore*/)
{
    tracker(getScore(dpValue), traceIter);
    tracker(dpValue);
}

template <typename TScoreValue, typename TTraceIter, typename TDPValue>
inline void track(Tracker<TScoreValue, TTraceIter, BandedChain> & tracker,
                  TDPValue const & dpValue,
                  TTraceIter const & /*traceIter*/,
                  False const & /*trackScore*/)
{
    tracker(dpValue);
}


template <typename TScoreValue, typename TTraceIter, typename TDpIter, typename TBand>
inline void
trackLastColumn(Tracker<TScoreValue, TTraceIter, BandedChain> & tracker,
                TDpIter & activeColIter,
                TDpIter const & beginDpColumn,
                TTraceIter & traceIter,
                TBand const & band,
                False const & /*switch*/)
{
        // only track last cell
    trackLastColumn(tracker, activeColIter, beginDpColumn, traceIter, band, True());
}

template <typename TScoreValue, typename TTraceIter, typename TDpIter>
inline void
trackLastColumn(Tracker<TScoreValue, TTraceIter, BandedChain> & tracker,
                TDpIter const & activeColIter,
                TDpIter const & beginDpColumn,
                TTraceIter & traceIter,
                Band<BandSwitchedOn<> > const & /*band*/,
                True const & /*switch*/)
{
    // need to scan last column from top to down to be conform with old overlap computation
    TDpIter it = beginDpColumn;
    traceIter -= (activeColIter - beginDpColumn);
    do
    {
        track(tracker, value(it), traceIter, True());
        ++it;
        ++traceIter;
    }while (it != activeColIter);
}

template <typename TScoreValue, typename TTraceIter, typename TDpIter>
inline void
trackLastColumn(Tracker<TScoreValue, TTraceIter, BandedChain> & tracker,
                TDpIter const & activeColIter,
                TDpIter const & beginDpColumn,
                TTraceIter & traceIter,
                Band<BandSwitchedOff> const & /*band*/,
                True const & /*switch*/)
{
    // need to scan last column from top to down to be conform with old overlap computation
    TDpIter it = beginDpColumn + tracker._posDim1;
    traceIter -= (activeColIter - (beginDpColumn + tracker._posDim1));
    do
    {
        track(tracker, value(it), traceIter, True());
        ++it;
        ++traceIter;
    }while (it != activeColIter);
}


// adapted initialization process of the banded chain alignment algorithm
template<typename TScoreValue, typename TTraceIterator, typename TDPIterator,
typename TSeqVIterator, typename TSeqHIterator, typename TDPFormula, typename TColumnType,
typename TAlignmentProfile>
inline void initializeMatrix(Tracker<TScoreValue, TTraceIterator, BandedChain> & tracker,
                             TTraceIterator & traceIter,
                             TDPIterator & activeColIter,
                             TSeqVIterator const & seqVBegin,
                             TSeqHIterator const & /*seqHIter*/,
                             TDPFormula const & /*dpFormula*/,
                             TColumnType const & /*columnType*/,
                             DPManager<TAlignmentProfile, Band<BandSwitchedOn> > & columnManager)
{
    typedef Tracker<TScoreValue, TTraceIterator, BandedChain> TTracker;
    typedef typename TTracker::TInitString TInitString;
    typedef typename Iterator<TInitString>::Type TInitStringIterator;

    TSeqVIterator seqVIter = seqVBegin + columnManager._spanSeqVBegin; // set begin of vertical seq relative to band
    TSeqVIterator seqVStop = seqVBegin + columnManager._spanSeqVEnd;  // set end of vertical sequence relative to band

    TInitStringIterator itInit = begin(tracker._initDim1Curr);
    // just copy the vertical init values of the tracker into the active column
    for (;seqVIter != seqVStop; ++seqVIter, ++itInit, ++traceIter, ++activeColIter)
    {
        value(activeColIter) = value(itInit);
        value(traceIter) = getOrigin(value(itInit));    // TODO(rmaerker): check how you can revise this
    }
    value(activeColIter) = value(itInit);
    value(traceIter) = getOrigin(value(itInit));
}

template<typename TScoreValue, typename TTraceIterator, typename TDPIterator,
typename TSeqVIterator, typename TSeqHIterator, typename TDPFormula, typename TColumnType,
typename TAlignmentProfile>
inline void initializeMatrix(Tracker<TScoreValue, TTraceIterator, BandedChain> & tracker,
                             TTraceIterator & traceIter,
                             TDPIterator & activeColIter,
                             TSeqVIterator const & seqVBegin,
                             TSeqHIterator const & /*seqHIter*/,
                             TDPFormula const & /*dpFormula*/,
                             TColumnType const & /*columnType*/,
                             DPManager<TAlignmentProfile, Band<BandSwitchedOff> > & columnManager)
{
    typedef Tracker<TScoreValue, TTraceIterator, BandedChain> TTracker;
    typedef typename TTracker::TInitString TInitString;
    typedef typename Iterator<TInitString>::Type TInitStringIterator;

    // we depend on the band size so we follow the seqVIterator
    TSeqVIterator seqVIter = seqVBegin + columnManager._spanSeqVBegin; // set begin of vertical seq relative to band
    TSeqVIterator seqVStop = seqVBegin + columnManager._spanSeqVEnd;  // set end of vertical sequence relative to band

    TInitStringIterator itInit = begin(tracker._initDim1Curr);
    value(activeColIter) = value(itInit);
    value(traceIter) = getOrigin(value(itInit));    // this needs to be revised some how.
    ++itInit;
    ++traceIter;
    ++activeColIter;
    // just copy the vertical init values of the tracker into the active column
    for( ;seqVIter != seqVStop; ++seqVIter, ++itInit, ++traceIter, ++activeColIter)
    {
        value(activeColIter) = value(itInit);
        value(traceIter) = getOrigin(value(itInit));
    }
    value(activeColIter) = value(itInit);
    value(traceIter) = getOrigin(value(itInit));
}

// TODO (rmaerker): plus one whenever we request the begin pos of a seed

template <typename TTraceSet, typename TInitSpotArray, typename TTracker, typename TSeed, typename TSeedSize, typename TGridPoint,
        typename TSeq0, typename TSeq1, typename TScoreScheme, typename TAlignmentProfile>
inline typename Value<TScoreScheme>::Type
_initializeBandedChain(TTraceSet & globalTraceSet,
        TInitSpotArray & initSpots,
        TTracker & tracker,
        TSeed const & seed,
        TSeedSize const & minBandWidth,
        TGridPoint & gridBegin,
        TGridPoint & gridEnd,
        TSeq0 const & seq0,
        TSeq1 const & seq1,
        TScoreScheme const & scoreSchemeGap,
        TScoreScheme const & scoreSchemeAnchor,
        TAlignmentProfile const & /*alignProfile*/)
{
    typedef typename Value<TInitSpotArray>::Type TInitSpot;
    typedef typename TInitSpot::T3 TDPValue;
    typedef typename Infix<TSeq0 const>::Type TInfix0;
    typedef typename Infix<TSeq1 const>::Type TInfix1;
    typedef typename Position<TSeed const>::Type TSeedPosition;
    typedef typename Value<TTraceSet>::Type TTraceArray;
    typedef typename Value<TScoreScheme>::Type TScore;

    // Infixes used for smaller grids
    TInfix0 infix0;
    TInfix1 infix1;

    TSeedPosition nextGridBeginDim0;
    TSeedPosition nextGridBeginDim1;

    TScore score = 0;

    // INITIALIZATION
    // TODO (rmaerker): change direction of dim0 and dim1 in seeds
    TSeedSize bandWidthDim0 = getBandShiftVertical(seed);
    TSeedSize bandWidthDim1 = getBandShiftHorizontal(seed);
    // define the values used for the initialization column and row
    appendValue(initSpots, TInitSpot(0, 0, TDPValue(0, +TraceBitMask::NONE)));
    for (unsigned i = 1; i <= _getBeginDim0(seed) + bandWidthDim0 + minBandWidth; ++i)
    {
        appendValue(initSpots, TInitSpot(i, 0, TDPValue(0, +TraceBitMask::NONE)));
    }
    for (unsigned j = 1; j <= _getBeginDim1(seed) + bandWidthDim1 + minBandWidth; ++j)
    {
        appendValue(initSpots, TInitSpot(0, j, TDPValue(0, +TraceBitMask::NONE)));
    }

    bool isCrossing = false;


    // if anchor crosses first row and column, then compute first anchor without preceeding standard dp grid
    if (isCrossingBeginDim0(_getBeginDim0(seed), minBandWidth) && isCrossingBeginDim1(_getBeginDim1(seed), minBandWidth))
    {
        isCrossing = true;
    }

    if (!isCrossing)
    {
        //####################################################
        //# part A: compute alignment to connect two anchors #
        //####################################################
        gridEnd.i1 = _getBeginDim0(seed) + minBandWidth + bandWidthDim0;
        gridEnd.i2 = _getBeginDim1(seed) + minBandWidth + bandWidthDim1;

        // define infixes covering the area which is to be computed
        infix0 = infix(seq0, gridBegin.i1, gridEnd.i1);
        infix1 = infix(seq1, gridBegin.i2, gridEnd.i2);

        // relative begin positions of next grid
        nextGridBeginDim0 = _max(0u, _getBeginDim0(seed) - minBandWidth);
        nextGridBeginDim1 = _max(0u, _getBeginDim1(seed) - minBandWidth);
        _init(tracker, gridEnd.i1 + 1, gridEnd.i2 + 1, nextGridBeginDim0, nextGridBeginDim1,
                ((minBandWidth << 1) + bandWidthDim0) + 1, ((minBandWidth << 1) + bandWidthDim1) + 1);
        _setInitSpots(tracker, initSpots);
        // compute the alignment
        score += _doComputeAlignment(globalTraceSet, tracker, initSpots, infix0, infix1,
                scoreSchemeGap, Band<BandSwitchedOff>(), TAlignmentProfile());
    }

    //###########################################
    //# part B: compute alignment around anchor #
    //###########################################

    // begin of next grid
    // TODO (rmaerker): use signed values here
    gridBegin.i1 = _max(0u, _getBeginDim0(seed) - minBandWidth);
    gridBegin.i2 = _max(0u, _getBeginDim1(seed) - minBandWidth);
    // end of next grid
    gridEnd.i1 = getEndDim0(seed) + minBandWidth;
    gridEnd.i2 = getEndDim1(seed) + minBandWidth;

    // define area covering infixes
    infix0 = infix(seq0, gridBegin.i1, gridEnd.i1);
    infix1 = infix(seq1, gridBegin.i2, gridEnd.i2);

    // reinitialize the tracker

    // begin of next grid to fill interspace after this anchor
    // relative begin positions of next grid
    nextGridBeginDim0 = getEndDim0(seed) - (minBandWidth + bandWidthDim0) - gridBegin.i1;
    nextGridBeginDim1 = getEndDim1(seed) - (minBandWidth + bandWidthDim1) - gridBegin.i2;

    // TODO (rmaerker): Caution with the length of tracker. What if it is not initialized
    if (!isCrossing)
    {
        _init(tracker, length(tracker._initDim0Next), length(tracker._initDim1Next), nextGridBeginDim0,
                nextGridBeginDim1, ((minBandWidth << 1) + bandWidthDim0) + 1, ((minBandWidth << 1) + bandWidthDim1) + 1);
    }
    else
    {
        _init(tracker, _getBeginDim0(seed) + minBandWidth + bandWidthDim0, _getBeginDim1(seed) + minBandWidth + bandWidthDim1, nextGridBeginDim0,
                nextGridBeginDim1, ((minBandWidth << 1) + bandWidthDim0) + 1, ((minBandWidth << 1) + bandWidthDim1) + 1);
    }

    _setInitSpots(tracker, initSpots);

    // compute the alignment of the anchor
    if (!isCrossing)
    {
        TTraceSet localTraceSet;
        score += _doComputeAlignment(
                localTraceSet, tracker, initSpots, infix0, infix1, scoreSchemeAnchor,
                Band<BandSwitchedOn<> >(((minBandWidth << 1) + bandWidthDim0), ((minBandWidth << 1) + bandWidthDim1)), TAlignmentProfile());
        _adaptLocalTracesToGlobalGrid(localTraceSet, gridBegin);
        glueTraces(globalTraceSet, localTraceSet);
    }
    else
    {
        score += _doComputeAlignment(
                globalTraceSet, tracker, initSpots, infix0, infix1, scoreSchemeAnchor,
                Band<BandSwitchedOn<> >(((minBandWidth << 1) + bandWidthDim0), ((minBandWidth << 1) + bandWidthDim1)),
                TAlignmentProfile());
    }
    return score;
}

template <typename TTraceSet, typename TInitSpotString, typename TTracker, typename TSeed, typename TSeedSize, typename TShiftSize, typename TGridPoint, typename TSeq0, typename TSeq1,
    typename TScoreScheme, typename TAlignmentProfile>
inline typename Value<TScoreScheme>::Type
_computeInterspaceArea(TTraceSet & globalTraces,
                    TInitSpotString & initSpots,
                    TTracker & tracker,
                    TSeed const & seed,
                    TSeedSize const & minBandWidth,
                    TShiftSize const & bandShiftHorizontal,
                    TShiftSize const & bandShiftVertical,
                    TGridPoint const & gridBegin,
                    TGridPoint const & gridEnd,
                    TSeq0 const & seq0,
                    TSeq1 const & seq1,
                    TScoreScheme const & scoreScheme,
                    TAlignmentProfile const & /*alignProfile*/)
{
    typedef typename Infix<TSeq0 const>::Type TInfix0;
    typedef typename Infix<TSeq1 const>::Type TInfix1;
    typedef typename Position<TSeq0>::Type TPos0;
    typedef typename Position<TSeq1>::Type TPos1;
    typedef Band<BandSwitchedOff> TBand;
    typedef typename Value<TScoreScheme>::Type TScoreValue;

    TTraceSet localTraces;
    // define area covering infixes
    TInfix0 infix0 = infix(seq0, gridBegin.i1, gridEnd.i1);
    TInfix1 infix1 = infix(seq1, gridBegin.i2, gridEnd.i2);

    // TODO (rmaerker): change direction of dim0 and dim1 in seeds

    // reinitialize the tracker
    TPos0 nextGridBeginDim0 = _getBeginDim0(seed) - minBandWidth - gridBegin.i1;
    TPos1 nextGridBeginDim1 = _getBeginDim1(seed) - minBandWidth - gridBegin.i2;

    _init(tracker, length(infix0) + 1, length(infix1) + 1, nextGridBeginDim0, nextGridBeginDim1,
            ((minBandWidth << 1) + bandShiftHorizontal) + 1, ((minBandWidth << 1) + bandShiftVertical) + 1);
    _setInitSpots(tracker, initSpots);

    // compute the alignment
    TScoreValue score = _doComputeAlignment(localTraces, tracker, initSpots, infix0, infix1, scoreScheme,
            TBand(), TAlignmentProfile());

    // in first matrix just copy the values of local traces to the global one
    _adaptLocalTracesToGlobalGrid(localTraces, gridBegin);
    glueTraces(globalTraces, localTraces);
    return score;
}

template <typename TTraceSet, typename TInitSpotString, typename TTracker, typename TGridPoint, typename TSeq0, typename TSeq1,
    typename TScoreScheme, typename TAlignmentProfile>
inline typename Value<TScoreScheme>::Type
_computeLastInterspaceArea(TTraceSet & globalTraces,
                    TInitSpotString & initSpots,
                    TTracker & tracker,
                    TGridPoint const & gridBegin,
                    TGridPoint const & gridEnd,
                    TSeq0 const & seq0,
                    TSeq1 const & seq1,
                    TScoreScheme const & scoreScheme,
                    TAlignmentProfile const & /*alignProfile*/)
{
    typedef typename Infix<TSeq0 const>::Type TInfix0;
    typedef typename Infix<TSeq1 const>::Type TInfix1;
    typedef typename Position<TSeq0>::Type TPos0;
    typedef typename Position<TSeq1>::Type TPos1;
    typedef typename Value<TScoreScheme>::Type TScoreValue;
    typedef String<typename TraceBitMask::Type> TTraceMatrix;

    TTraceSet localTraceSet;
    // define area covering infixes
    TInfix0 infix0 = infix(seq0, gridBegin.i1, gridEnd.i1);
    TInfix1 infix1 = infix(seq1, gridBegin.i2, gridEnd.i2);


    _init(tracker, length(infix0) + 1, length(infix1) + 1, 0, 0, length(infix0) + 1, length(infix1) + 1);
    _setInitSpots(tracker, initSpots);

    // compute the alignment
    Band<BandSwitchedOff> band;

    TTraceMatrix traceMatrix;
    resize(traceMatrix, getTraceMatrixSize(infix0, infix1, band, True()), +TraceBitMask::NONE);

    computeDPMatrix(tracker, traceMatrix, infix0, infix1, scoreScheme, band, TAlignmentProfile());

    // DEBUG CODE
    // TODO (rmaerker): implement this
//    std::cout << seq0 << std::endl;
//    std::cout << seq1 << std::endl;
//  std::cout << "Traceback:\t";
//  for (unsigned i = 0; i < length(tracker._tracePoints);++i)
//  {
//      std::cout << tracker._tracePoints[i] - static_cast<TTraceIterator>(begin(traceMatrix)) << std::endl;
//      std::cout << _determineSeq0Pos(tracker._tracePoints[i] - static_cast<TTraceIterator>(begin(traceMatrix)), seq0, seq1, band) << "\t";
//      std::cout << _determineSeq1Pos(tracker._tracePoints[i] - static_cast<TTraceIterator>(begin(traceMatrix)), seq0, seq1, band) << "\t";
//      std::cout << "Trace Value: " << (int) value(tracker._tracePoints[i]) << std::endl;
//
//  }

    computeTraceback(localTraceSet, tracker, traceMatrix, infix0, infix1, band, TAlignmentProfile());

    // in first matrix just copy the values of local traces to the global one
    _adaptLocalTracesToGlobalGrid(localTraceSet, gridBegin);
    glueTraces(globalTraces, localTraceSet);

    // TODO (rmaerker): check the banded trace back version. does not work correctly.
    // DEBUG CODE

//  std::cout << "Init Dim0:\t";
//  for (unsigned i = 0; i < length(tracker._initDim0Next);++i)
//  {
//      std::cout << tracker._initDim0Next[i] << "\t";
//  }
//  std::cout << "\nInit Dim1:\t";
//  for (unsigned i = 0; i < length(tracker._initDim1Next);++i)
//  {
//      std::cout << tracker._initDim1Next[i] << "\t";
//  }
//  std::cout << "\nInit Spots:\t";
//  for (unsigned i = 0; i < length(initSpots);++i)
//  {
//      std::cout << initSpots[i] << "\t";
//  }
//  std::cout << "\n";
    return tracker._maxScore;
}

template <typename TTraceSet, typename TInitSpotString, typename TTracker, typename TSeed, typename TSeedSize, typename TShiftSize, typename TGridPoint, typename TSeq0, typename TSeq1,
    typename TScoreScheme, typename TAlignmentProfile>
inline typename Value<TScoreScheme>::Type
_computeAnchorArea(TTraceSet & globalTraces,
                   TInitSpotString & initSpots,
                   TTracker & tracker,
                   TSeed const & seed,
                   TSeedSize const & minBandWidth,
                   TShiftSize const & bandShiftHorizontal,
                   TShiftSize const & bandShiftVertical,
                   TGridPoint const & gridBegin,
                   TGridPoint const & gridEnd,
                   TSeq0 const & seq0,
                   TSeq1 const & seq1,
                   TScoreScheme const & scoreScheme,
                   TAlignmentProfile const & alignProfile)
{
    typedef typename Infix<TSeq0 const>::Type TInfix0;
    typedef typename Infix<TSeq0 const>::Type TInfix1;
    typedef typename Position<TSeq0 const>::Type TPos0;
    typedef typename Position<TSeq1 const>::Type TPos1;
    typedef Band<BandSwitchedOn<> > TBand;
    typedef typename Value<TScoreScheme>::Type TScore;

    TTraceSet localTraces;
    // define area covering infixes
    TInfix0 infix0 = infix(seq0, gridBegin.i1, gridEnd.i1);
    TInfix1 infix1 = infix(seq1, gridBegin.i2, gridEnd.i2);

    // TODO (rmaerker): change direction of dim0 and dim1 in seeds

    // reinitialize the tracker
    TPos0 nextGridBeginDim0 = getEndDim0(seed) - (minBandWidth + bandShiftHorizontal) - gridBegin.i1;
    TPos1 nextGridBeginDim1 = getEndDim1(seed) - (minBandWidth + bandShiftVertical) - gridBegin.i2;

    _init(tracker, length(tracker._initDim0Next),length(tracker._initDim1Next), nextGridBeginDim0, nextGridBeginDim1,
            ((minBandWidth << 1) + bandShiftHorizontal) + 1, ((minBandWidth << 1) + bandShiftVertical) + 1);
    _setInitSpots(tracker, initSpots);

    // compute the alignment
    TScore score = _doComputeAlignment(localTraces, tracker, initSpots, infix0, infix1, scoreScheme,
            TBand(((minBandWidth << 1) + bandShiftHorizontal), ((minBandWidth << 1) + bandShiftVertical)), alignProfile);

    // in first matrix just copy the values of local traces to the global one
    _adaptLocalTracesToGlobalGrid(localTraces, gridBegin);
    glueTraces(globalTraces, localTraces);
    return score;
}

// The LAGAN algorithms defines a rough global map of the alignment by identifying high scoring anchors.
// In the followed banded chain alignment two anchors are connect by filling the gap between them using
// a standard dp algorithm.
template <typename TTraceSet, typename TTracker, typename TInitSpots, typename TSequence0, typename TSequence1, typename TScoreScheme,
        typename TBand>
inline typename Value<TScoreScheme>::Type
_doComputeAlignment(TTraceSet & localTraceSet,
                   TTracker & tracker,
                   TInitSpots & initSpots,
                   TSequence0 const & seq0,
                   TSequence1 const & seq1,
                   TScoreScheme const & scoreScheme,
                   TBand const & band,
                   AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn> const & profile)
{
    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn> TAlignmentProfile;
    typedef typename Value<TScoreScheme>::Type TScoreValue;
    typedef String<typename TraceBitMask::Type> TTraceMatrix;
    typedef typename Iterator<TTraceMatrix>::Type TTraceIterator;
    typedef typename Value<TScoreScheme>::Type TScore;

    // handle empty sequences
    if (empty(seq0) || empty(seq1) )
    {
        SEQAN_ASSERT_FAIL("Trying to compute alignment with empty sequences in banded chain alignment!");
    }

    TTraceMatrix traceMatrix;
    resize(traceMatrix, getTraceMatrixSize(seq0, seq1, band, True()), +TraceBitMask::NONE);

    computeDPMatrix(tracker, traceMatrix, seq0, seq1, scoreScheme, band, profile);

    // DEBUG CODE
    // TODO (rmaerker): implement this
    std::cout << seq0 << std::endl;
    std::cout << seq1 << std::endl;
    std::cout << "Traceback:\t";
    for (unsigned i = 0; i < length(tracker._tracePoints);++i)
    {
        std::cout << tracker._tracePoints[i] - static_cast<TTraceIterator>(begin(traceMatrix)) << std::endl;
        std::cout << _determineSeq0Pos(tracker._tracePoints[i] - static_cast<TTraceIterator>(begin(traceMatrix)), seq0, seq1, band) << "\t";
        std::cout << _determineSeq1Pos(tracker._tracePoints[i] - static_cast<TTraceIterator>(begin(traceMatrix)), seq0, seq1, band) << "\t";
        std::cout << "Trace Value: " << (int) value(tracker._tracePoints[i]) << std::endl;

    }

    computeTraceback(localTraceSet, initSpots, tracker, traceMatrix, seq0, seq1, band, profile);

    // TODO (rmaerker): check the banded trace back version. does not work correctly.
    // DEBUG CODE

    std::cout << "Init Dim0:\t";
    for (unsigned i = 0; i < length(tracker._initDim0Next);++i)
    {
        std::cout << tracker._initDim0Next[i] << "\t";
    }
    std::cout << "\nInit Dim1:\t";
    for (unsigned i = 0; i < length(tracker._initDim1Next);++i)
    {
        std::cout << tracker._initDim1Next[i] << "\t";
    }
    std::cout << "\nInit Spots:\t";
    for (unsigned i = 0; i < length(initSpots);++i)
    {
        std::cout << initSpots[i] << "\t";
    }
    std::cout << "\n";
    return tracker._maxScore;
}


// The LAGAN algorithms defines a rough global map of the alignment by identifying high scoring anchors.
// In the followed banded chain alignment two anchors are connect by filling the gap between them using
// a standard dp algorithm.
template <typename TTraceSet, typename TSeedSet, typename TSequence0, typename TSequence1, typename TScoreScheme, typename TBandSize>
inline typename Value<TScoreScheme>::Type
_computeAlignment(TTraceSet & globalTraceSet,
                TSeedSet const & seedSet,
                TSequence0 const & seq0,
                TSequence1 const & seq1,
                TScoreScheme const & scoreSchemeGap,
                TScoreScheme const & scoreSchemeAnchor,
                TBandSize const & minBandWidth,
                AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn> const & profile)
{
    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn> TAlignmentProfile;

    typedef typename Position<TSequence0>::Type TPos0;
    typedef typename Position<TSequence1>::Type TPos1;
    typedef Pair<TPos0, TPos1> TGridPoint;

    typedef typename Value<TScoreScheme>::Type TScoreValue;
    typedef DPValue<TScoreValue, AffineGaps> TDPValue;

    typedef typename Iterator<TSeedSet const>::Type TSeedSetIterator;
    typedef typename Value<TSeedSet>::Type TSeed;
    typedef typename Position<TSeed>::Type TSeedPosition;
    typedef typename Size<TSeed>::Type TSeedSize;
    typedef TraceSegment<TSeedPosition, TSeedSize> TTraceSegment;
    typedef String<TTraceSegment> TTraceSegments;

    typedef typename Iterator<String<typename TraceBitMask::Type>, Standard>::Type TTraceIterator;
    typedef typename Iterator<TTraceSegments>::Type TTraceSegmentsIterator;
    typedef Tracker<TScoreValue, TTraceIterator, BandedChain> TTracker;
    typedef Triple<TPos0, TPos1, TDPValue> TInitSpot;
    typedef String<TInitSpot> TInitSpots;

    // handle cases of empty sequences
    SEQAN_ASSERT_MSG(!empty(seq0), "Cannot compute alignment on empty sequence (horizontal sequence)");
    SEQAN_ASSERT_MSG(!empty(seq1), "Cannot compute alignment on empty sequence (vertical sequence)");

    // handle case of empty seed set
    if (length(seedSet) < 1)
    {
        resize(globalTraceSet, 1, Exact());
        appendValue(globalTraceSet[0], TTraceSegment(length(seq0),0, length(seq1), +TraceBitMask::VERTICAL));
        appendValue(globalTraceSet[0], TTraceSegment(0,0, length(seq0), +TraceBitMask::HORIZONTAL));
        return 0;
    }

    // PREPROCESSING

    // initialize grid points


    // initiailze seedSet iterator
    TSeedSetIterator it = findFirstAnchor(seedSet, minBandWidth);
    TSeedSetIterator itEnd = findLastAnchor(seedSet, seq0, seq1, minBandWidth);

    SEQAN_ASSERT(itEnd != static_cast<TSeedSetIterator>(end(seedSet)));

    // handle the first anchor
    TSeed  seed = value(it);

    TGridPoint gridBegin(0,0);
    TGridPoint gridEnd(0,0);

    // define array to keep track of the initialization spots for the next dp grid
    TInitSpots initSpots;

    TTracker tracker;

    TScoreValue tmp;
    TScoreValue score = 0;
    // INITIALIZATION
    tmp = _initializeBandedChain(globalTraceSet, initSpots, tracker, seed, minBandWidth,
            gridBegin, gridEnd, seq0, seq1, scoreSchemeGap, scoreSchemeAnchor, profile);
    score += tmp;


    // MAIN
    TSeedSize bandShiftHorizontal = getBandShiftVertical(seed);
    TSeedSize bandShiftVertical = getBandShiftHorizontal(seed);
    while (it != itEnd)
    {
        //####################################################
        //# part A: compute alignment to connect two anchors #
        //####################################################
        // begin point of grid
        gridBegin.i1  = getEndDim0(seed) - minBandWidth - bandShiftHorizontal; // - band shift
        gridBegin.i2  = getEndDim1(seed) - minBandWidth - bandShiftVertical; // - band shift

        // move to next anchor
        seed = value(++it);

        //new band shift
        bandShiftHorizontal = getBandShiftVertical(seed);
        bandShiftVertical = getBandShiftHorizontal(seed);

        // end point of grid
        gridEnd.i1 = _getBeginDim0(seed) + minBandWidth + bandShiftHorizontal; // + bandShift
        gridEnd.i2 = _getBeginDim1(seed) + minBandWidth + bandShiftVertical; // + bandShift

        tmp = _computeInterspaceArea(globalTraceSet, initSpots, tracker,
                seed, minBandWidth, bandShiftHorizontal, bandShiftVertical, gridBegin, gridEnd, seq0, seq1, scoreSchemeGap, profile);
        score += tmp;

        //###########################################
        //# part B: compute alignment around anchor #
        //###########################################

        // define the grid for the anchor
        gridBegin.i1  = _getBeginDim0(seed) - minBandWidth;
        gridBegin.i2  = _getBeginDim1(seed) - minBandWidth;
        gridEnd.i1 = getEndDim0(seed) + minBandWidth;  // carfully don't let it exceed the border
        gridEnd.i2 = getEndDim1(seed) + minBandWidth;

        tmp = _computeAnchorArea(globalTraceSet, initSpots, tracker,
                seed, minBandWidth, bandShiftHorizontal, bandShiftVertical, gridBegin, gridEnd, seq0, seq1, scoreSchemeAnchor, profile);
        score += tmp;
    }

    // POSTPROCESSING
    // handle last anchor - can be crossing last row and/or column

    //####################################################
    //# part A: compute alignment to connect two anchors #
    //####################################################

    // begin point of grid
    gridBegin.i1  = getEndDim0(seed) - minBandWidth - bandShiftHorizontal;    // ensure negative values are allowed
    gridBegin.i2  = getEndDim1(seed) - minBandWidth - bandShiftVertical;

    // move to next anchor
    seed = value(++it);

    //new band shift
    bandShiftHorizontal = getBandShiftVertical(seed);
    bandShiftVertical = getBandShiftHorizontal(seed);

    // end point of grid
    gridEnd.i1 = _getBeginDim0(seed) + minBandWidth + bandShiftHorizontal;
    gridEnd.i2 = _getBeginDim1(seed) + minBandWidth + bandShiftVertical;

    tmp = _computeInterspaceArea(globalTraceSet, initSpots, tracker,
            seed, minBandWidth, bandShiftHorizontal, bandShiftVertical, gridBegin, gridEnd, seq0, seq1, scoreSchemeGap, TAlignmentProfile());
    score += tmp;

    //###########################################
    //# part B: compute alignment around anchor #
    //###########################################


    // define the grid for the anchor
    gridBegin.i1  = _getBeginDim0(seed) - minBandWidth;
    gridBegin.i2  = _getBeginDim1(seed) - minBandWidth;
    bool isCrossing = false;
    if (isCrossingEndDim0(getEndDim0(seed), minBandWidth, seq0) && isCrossingEndDim1(getEndDim1(seed), minBandWidth, seq1))
    {
        isCrossing = true;

    }
    gridEnd.i1 = _min(length(seq0), getEndDim0(seed) + minBandWidth);  // carefully don't let it exceed the border
    gridEnd.i2 = _min(length(seq1), getEndDim1(seed) + minBandWidth);

    tmp = _computeAnchorArea(globalTraceSet, initSpots, tracker,
            seed, minBandWidth, bandShiftHorizontal, bandShiftVertical, gridBegin, gridEnd, seq0, seq1, scoreSchemeAnchor, TAlignmentProfile());
    score += tmp;

    if (!isCrossing)
    { // compute last rectangle
        // begin point of grid
        gridBegin.i1  = getEndDim0(seed) - minBandWidth - bandShiftHorizontal;    // ensure negative values are allowed
        gridBegin.i2  = getEndDim1(seed) - minBandWidth - bandShiftVertical;
        // end point of grid
        gridEnd.i1 = length(seq0);
        gridEnd.i2 = length(seq1);

        // compute the alignment
        _computeLastInterspaceArea(globalTraceSet, initSpots, tracker,
                    gridBegin, gridEnd, seq0, seq1, scoreSchemeGap, TAlignmentProfile());
    }
    std::cout << "The summed score: " << score << std::endl;
    return tracker._maxScore;
}

// use 15 as default band size such as is in the LAGAN paper
template <typename TPosition, typename TSize, typename TSequenceH, typename TSequenceV, typename TSeedSpec,
        typename TSeedSetSpec, typename TSeedConfig, typename TScoreScheme, typename TBandSize>
typename Value<TScoreScheme>::Type
_runBandedChainAlignment(StringSet<String<TraceSegment<TPosition, TSize> > > & globalTraceSet,
                         TSequenceH const & seq0,
                         TSequenceV const & seq1,
                         SeedSet<TSeedSpec, TSeedSetSpec, TSeedConfig> const & seedSet,
                         TScoreScheme const & scoreSchemeGap,
                         TScoreScheme const & scoreSchemeAnchor,
                         TBandSize const & minBandWidth)
{
    typedef AlignmentProfile<Global<BandedChain>, AffineGaps, TracebackSwitchedOn > TAlignmentProfile;
    typedef TraceSegment<TPosition, TSize> TTraceSegment;
    typedef String<TTraceSegment> TTraceSegArray;
    typedef StringSet<TTraceSegArray> TTraceSegArraySet;
    typedef typename Iterator<TTraceSegArray, Standard>::Type TTraceSegArrayIterator;
    typedef typename Value<TScoreScheme>::Type TScoreValue;

    TScoreValue score = _computeAlignment(globalTraceSet, seedSet, seq0, seq1, scoreSchemeGap, scoreSchemeAnchor,
                                          minBandWidth, TAlignmentProfile());

    // record overlaps at beginning or/and at end of sequences if any
    for (unsigned i = 0; i < length(globalTraceSet); ++i)
    {
        // #1 check if the first trace is overlapping
        TTraceSegArrayIterator iter = end(globalTraceSet[i]) - 1;
        if (getBeginHorizontal(value(iter)) != 0)
        { // fill horizontal gaps at beginning
            SEQAN_ASSERT_EQ_MSG(getBeginVertical(value(iter)), 0u,
                    "Invalid trace segment at beginning of matrix (last recorded trace segment)!");
            recordSegment(globalTraceSet[i], 0, 0, getBeginHorizontal(value(iter)), +TraceBitMask::HORIZONTAL);
        }
        if (getBeginVertical(value(iter)) != 0)
        { // fill vertical gaps at beginning
            SEQAN_ASSERT_EQ_MSG(getBeginHorizontal(value(iter)), 0u,
                    "Invalid trace segment at beginning of matrix (last recorded trace segment)!");
            recordSegment(globalTraceSet[i], 0, 0, getBeginVertical(value(iter)), +TraceBitMask::VERTICAL);
        }
        // #2 check if the last trace is overlapping
        iter = begin(globalTraceSet[i]);
        if (getEndHorizontal(value(iter)) != length(seq0))
        { // fill horizontal gaps at end
            SEQAN_ASSERT_EQ_MSG(getEndVertical(value(iter)), length(seq1),
                    "Invalid trace segment at end of matrix (first recorded trace segment)!");
            // insert at the beginning ...
            insertValue(globalTraceSet[i], 0, TTraceSegment(getEndHorizontal(value(iter)), getEndVertical(value(iter)),
                    length(seq0) - getEndHorizontal(value(iter)), +TraceBitMask::HORIZONTAL));
        }
        if (getEndVertical(value(iter)) != length(seq1))
        { // fill vertical gaps at end
            SEQAN_ASSERT_EQ_MSG(getEndHorizontal(value(iter)), length(seq0),
                    "Invalid trace segment at end of matrix (first recorded trace segment)!");
            // insert at the beginning ...
            insertValue(globalTraceSet[i], 0, TTraceSegment(getEndHorizontal(value(iter)), getEndVertical(value(iter)),
                    length(seq1) - getEndVertical(value(iter)), +TraceBitMask::VERTICAL));
        }
    }
    return score;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_SEEDS2_BANDED_CHAIN_ALIGNMENT_IMPL_H_
