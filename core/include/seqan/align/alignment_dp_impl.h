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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_IMPL_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_IMPL_H_

namespace seqan
{

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

// ----------------------------------------------------------------------------
// Function _checkBandProperties()
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TSequenceV, typename TAlignmentProfile>
inline bool _checkBandProperties(TSequenceH const & /*seqH*/,
                                 TSequenceV const & /*seqV*/,
                                 Band<BandSwitchedOff> const & /*band*/,
                                 TAlignmentProfile const & /*alignProfile*/)
{
    return true;
}

template <typename TSequenceH, typename TSequenceV, typename TBandSpec, typename TAlignmentProfile>
inline bool _checkBandProperties(TSequenceH const & seqH,
                                 TSequenceV const & seqV,
                                 Band<BandSwitchedOn<TBandSpec> > const & band,
                                 TAlignmentProfile const & /*alignProfile*/)
{
    typedef typename MakeSigned< typename Size<TSequenceH>::Type >::Type TSignedSize;

    // Check if the intersection between band and DP matrix is empty.
    if (getUpperDiagonal(band) < (0 - static_cast<TSignedSize>(length(seqV)))
            || getLowerDiagonal(band) > static_cast<TSignedSize>(length(seqH)))
    {
        return false;
    }

    // If the band begins before the beginning of the horizontal sequence
    // then check if free end-gaps are enabled at the beginning of the vertical sequence.
    if (getUpperDiagonal(band) < 0 && !IsFreeEndGap<TAlignmentProfile, FirstColumn>::VALUE)
        return false;

    // If the band begins before the beginning of the vertical sequence
    // then check if free end-gaps are enabled at the beginning of the horizontal sequence.
    if (getLowerDiagonal(band) > 0 && !IsFreeEndGap<TAlignmentProfile, FirstRow>::VALUE)
        return false;

    // If the band ends behind the end of the vertical sequence
    // then check if free end-gaps are enabled at the end of the horizontal sequence.
    if (getUpperDiagonal(band) + static_cast<TSignedSize>(length(seqV)) < static_cast<TSignedSize>(length(seqH))
            && !IsFreeEndGap<TAlignmentProfile, LastRow>::VALUE)
    {
        return false;
    }

    // If the band ends behind the end of the horizontal sequence
    // then check if free end-gaps are enabled at the end of the vertical sequence.
    if (getLowerDiagonal(band) + static_cast<TSignedSize>(length(seqV)) > static_cast<TSignedSize>(length(seqH))
            && !IsFreeEndGap<TAlignmentProfile, LastColumn>::VALUE)
    {
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------
// Function _invalidDPSettings()
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TSequenceV, typename TBand, typename TAlignmentProfile>
inline bool _isValidDPSettings(TSequenceH const & seqH,
                               TSequenceV const & seqV,
                               TBand const & band,
                               TAlignmentProfile const & alignProfile)
{
    // Check if the sequences are empty.
    if (empty(seqH) || empty(seqV))
    {
        return false;
    }

    return _checkBandProperties(seqH, seqV, band, alignProfile);
}

// ----------------------------------------------------------------------------
// function _setInitValue()
// ----------------------------------------------------------------------------

template <typename TDPIter, typename TTracker, typename TBand, typename TColumnType>
inline void
_setInitValue(TDPIter & /*activeColIter*/,
               TTracker const & /*tracker*/,
               TBand const & /*band*/,
               TColumnType const /*columnType*/)
{
    // nothing to do
}

// ----------------------------------------------------------------------------
// function _getHorizontalValue()
// ----------------------------------------------------------------------------

template <typename TDPValue, typename TDPIter, typename TColumnType>
inline void
_getHorizontalValue(TDPValue & /*prevH*/, TDPIter const & /*prevIter*/, TColumnType const & /*columnType*/)
{
    // nothing to do
}

template <typename TDPValue, typename TDPIter>
inline void
_getHorizontalValue(TDPValue & prevH, TDPIter & prevIter, DPWideBandLastPhase const & /*columnType*/)
{
   prevH = value(++prevIter);
}

template <typename TDPValue, typename TDPIter>
inline void
_getHorizontalValue(TDPValue & prevH, TDPIter & prevIter, DPBandLastPhase const & /*columnType*/)
{
   prevH = value(++prevIter);
}

// ----------------------------------------------------------------------------
// Function _beginBandFirstPhase()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline typename Iterator<TSequenceH const>::Type
_beginBandFirstPhase(TSequenceH  const & seqH,
                     TSequenceV const & /*seqV*/,
                     Band<TBandSpec> const & band)
{
    typedef typename Size<Band<TBandSpec> >::Type TSize;
    return begin(seqH) + _min(static_cast<TSize>(length(seqH)), _max(0, getLowerDiagonal(band)));
}

// ----------------------------------------------------------------------------
// Function _beginWideBandFirstPhase()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline typename Iterator<TSequenceH const>::Type
_beginWideBandFirstPhase(TSequenceH const & seqH,
                         TSequenceV const & /*seqV*/,
                         Band<TBandSpec> const & band)
{
    typedef typename Size<Band<TBandSpec> >::Type TSize;
    return begin(seqH) + _min(static_cast<TSize>(length(seqH)), _max(0, getLowerDiagonal(band)));
}

// ----------------------------------------------------------------------------
// Function _beginBandMiddlePhase()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline typename Iterator<TSequenceH const>::Type
_beginBandMiddlePhase(TSequenceH const & seqH,
                  TSequenceV const & /*seqV*/,
                  Band<TBandSpec> const & band)
{
    typedef typename Size<Band<TBandSpec> >::Type TSize;
    return begin(seqH) + _min(static_cast<TSize>(length(seqH)), _max(0, getUpperDiagonal(band)));
}

// ----------------------------------------------------------------------------
// Function _beginWideBandMiddlePhase()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline typename Iterator<TSequenceH const>::Type
_beginWideBandMiddlePhase(TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Band<TBandSpec> const & band)
{
    typedef typename Size<Band<TBandSpec> >::Type TSize;
    return begin(seqH) + _min(static_cast<TSize>(length(seqH)),
                              _max(0, (static_cast<TSize>(length(seqV)) + getLowerDiagonal(band))));
}

// ----------------------------------------------------------------------------
// Function _beginBandEndPhase()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline typename Iterator<TSequenceH const>::Type
_beginBandEndPhase(TSequenceH  const & seqH,
                 TSequenceV const & seqV,
                 Band<TBandSpec> const & band)
{
    typedef typename Size<Band<TBandSpec> >::Type TSize;
    return begin(seqH) + _min(static_cast<TSize>(length(seqH)),
                              static_cast<TSize>(length(seqV)) + getLowerDiagonal(band));
}

// ----------------------------------------------------------------------------
// Function _beginWideBandEndPhase()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline typename Iterator<TSequenceH const>::Type
_beginWideBandEndPhase(TSequenceH const & seqH,
                     TSequenceV const & /*seqV*/,
                     Band<TBandSpec> const & band)
{
    typedef typename Size<Band<TBandSpec> >::Type TSize;
    return begin(seqH) + _min(static_cast<TSize>(length(seqH)), _max(0, getUpperDiagonal(band)));
}

// ----------------------------------------------------------------------------
// Function _endBand()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline typename Iterator<TSequenceH const>::Type
_endBand(TSequenceH const & seqH,
         TSequenceV const & seqV,
         Band<TBandSpec> const & band)
{
    typedef typename Size<Band<TBandSpec> >::Type TSize;
    return  begin(seqH) + _min(static_cast<TSize>(length(seqH)), getUpperDiagonal(band) + static_cast<TSize>(length(seqV)));
}


// ----------------------------------------------------------------------------
// Function _endWideBand()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline typename Iterator<TSequenceH const>::Type
_endWideBand(TSequenceH const & seqH,
             TSequenceV const & seqV,
             Band<TBandSpec> const & band)
{
    typedef typename Size<Band<TBandSpec> >::Type TSize;
    return  begin(seqH) + _min(static_cast<TSize>(length(seqH)), getUpperDiagonal(band) + static_cast<TSize>(length(seqV)));
}

// ----------------------------------------------------------------------------
// function increment() - conitional increment
// ----------------------------------------------------------------------------

template <typename TIterator>
inline void increment(TIterator & iter, True const & /*traceSpec*/)
{
    ++iter;
}

template <typename TIterator>
inline void increment(TIterator & /*iter*/, False const & /*traceSpec*/)
{
        // nothing to do
}

// ----------------------------------------------------------------------------
// function track()
// ----------------------------------------------------------------------------

template <typename TTracker, typename TDPValue, typename TTraceIterator>
inline void track(TTracker & tracker,
                  TDPValue const & dpValue,
                  TTraceIterator const & traceIter,
                  True const & /*trackScore*/)
{
    tracker(getScore(dpValue), traceIter);
}

template <typename TTracker, typename TDPValue, typename TTraceIterator>
inline void track(TTracker & /*traceTracker*/,
                  TDPValue const & /*dpValue*/,
                  TTraceIterator const & /*traceIter*/,
                  False const & /*trackScore*/)
{
    //nothing to do per default
}

// ----------------------------------------------------------------------------
// function trackLastColumn()
// ----------------------------------------------------------------------------

template <typename TTracker, typename TDpIter, typename TTraceIter, typename TBand>
inline void
trackLastColumn(TTracker & tracker,
                TDpIter & activeColIter,
                TDpIter const & /*beginDpColumn*/,
                TTraceIter & traceIter,
                TBand const & /*band*/,
                False const & /*switch*/)
{
    // only track last cell
    track(tracker, value(activeColIter), traceIter, True());
}

template <typename TTracker, typename TDpIter, typename TTraceIter, typename TBand>
inline void
trackLastColumn(TTracker & tracker,
                 TDpIter const & activeColIter,
                 TDpIter const & beginDpColumn,
                 TTraceIter & traceIter,
                 TBand const & /*band*/,
                 True const & /*switch*/)
{
    // need to scan last column from top to down to be conform with old overlap computation
    TDpIter it = beginDpColumn;
    traceIter -= activeColIter - beginDpColumn; // set traceback iterator to begin of column
    while (it != activeColIter)
    {
        track(tracker, value(it), traceIter, True());
        ++it;
        ++traceIter;
    }
    track(tracker, value(activeColIter), traceIter, True()); // track last cell
}

// ----------------------------------------------------------------------------
// function computeCell()
// ----------------------------------------------------------------------------

template<typename TTracker, typename TTraceIterator, typename TDPValue, typename TSeq0Value,
        typename TSeq1Value, typename TRecursionFormula, typename TCell>
inline void computeCell(TTracker & tracker,
        TTraceIterator & traceIter,
        TDPValue & activeCell,
        TDPValue const & previousDiagonal,
        TDPValue const & previousHorizontal,
        TDPValue const & previousVertical,
        TSeq0Value const & seqHValue,
        TSeq1Value const & seqVValue,
        TRecursionFormula const & recursionFormula,
        TCell const & /*cellSpec*/,
        True const & /*traceSpec*/)
{
    //call functor to compute the score of the active cell and store the resulting trace
    value(traceIter) = recursionFormula(activeCell, previousDiagonal, previousHorizontal, previousVertical, seqHValue,
            seqVValue, typename GetDpDirection<TCell>::Type());
    track(tracker, activeCell, traceIter, typename IsToTrack<TCell>::Type());
}

template<typename TTracker, typename TTraceIterator, typename TDPValue, typename TSeq0Value,
        typename TSeq1Value, typename TRecursionFormula, typename TCell>
inline void computeCell(TTracker & tracker,
        TTraceIterator & traceIter,
        TDPValue & activeCell,
        TDPValue const & previousDiagonal,
        TDPValue const & previousHorizontal,
        TDPValue const & previousVertical,
        TSeq0Value const & seqHValue,
        TSeq1Value const & seqVValue,
        TRecursionFormula const & recursionFormula,
        TCell const & /*cellSpec*/,
        False const & /*traceSpec*/)
{
    //call functor to compute the score of the active cell and store the resulting trace
    recursionFormula(activeCell, previousDiagonal, previousHorizontal, previousVertical, seqHValue, seqVValue,
                     typename GetDpDirection<TCell>::Type());
    track(tracker, activeCell, traceIter, typename IsToTrack<TCell>::Type());
}

//// ----------------------------------------------------------------------------
//// Function _computeFirstCell()
//// ----------------------------------------------------------------------------
//
//template <typename TTraceTracker, typename TTraceIterator, typename TDPIterator, typename TDPValue, typename TSeqVIterator,
//    typename TSeqHValue, typename TDPFormula, typename TAlignProfile, typename TBandSpec>
//inline bool _computeFirstCell(TTraceTracker & tracker,
//                              TTraceIterator & traceIter,
//                              TDPIterator & dpIter,
//                              TDPValue const & prevD,
//                              TDPValue const & prevH,
//                              TDPValue const & prevV,
//                              TSeqVIterator const & seqVIter,
//                              TSeqVIterator const & seqVStop,
//                              TSeqHValue const & seqHValue,
//                              TDPFormula const & dpFormula,
//                              DPBandClosurePhase const & /*dpPhase*/,
//                              DPManager<TAlignProfile, Band<TBandSpec> > const & /*dpManager*/)
//{
//    typedef typename SetUpColumnManager<TAlignProfile, DPBandClosurePhase>::Type TColumnManager;
//    typedef Band<TBandSpec> TBand;
//    typedef typename IsTracebackOn<TAlignProfile>::Type TIsTracebackOn;
//
//
//    if (seqVIter == seqVStop)  // if true than first cell is also last cell and must be tracked.
//    {
//        computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
//                    Cell<typename GetDpDirection<typename TColumnManager::TFirstCell>::Type, True>(), TIsTracebackOn());
//        return true;
//    }
//    computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
//                 typename TColumnManager::TFirstCell(), TIsTracebackOn());
//    return false;
//}
//
//template <typename TTraceTracker, typename TTraceIterator, typename TDPIterator, typename TDPValue, typename TSeqVIterator,
//    typename TSeqHValue, typename TDPFormula, typename TAlignProfile, typename TBandSpec>
//inline bool _computeFirstCell(TTraceTracker & tracker,
//                              TTraceIterator & traceIter,
//                              TDPIterator & dpIter,
//                              TDPValue const & prevD,
//                              TDPValue const & prevH,
//                              TDPValue const & prevV,
//                              TSeqVIterator const & seqVIter,
//                              TSeqVIterator const & seqVStop,
//                              TSeqHValue const & seqHValue,
//                              TDPFormula const & dpFormula,
//                              DPWideBandClosurePhase const & /*dpPhase*/,
//                              DPManager<TAlignProfile, Band<TBandSpec> > const & /*dpManager*/)
//{
//    typedef typename SetUpColumnManager<TAlignProfile, DPWideBandClosurePhase>::Type TColumnManager;
//    typedef Band<TBandSpec> TBand;
//    typedef typename IsTracebackOn<TAlignProfile>::Type TIsTracebackOn;
//
//
//    if (seqVIter == seqVStop)  // if true than first cell is also last cell and must be tracked.
//    {
//        computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
//                    Cell<typename GetDpDirection<typename TColumnManager::TFirstCell>::Type, True>(), TIsTracebackOn());
//        return true;
//    }
//    computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
//                 typename TColumnManager::TFirstCell(), TIsTracebackOn());
//    return false;
//}
//
//template <typename TTraceTracker, typename TTraceIterator, typename TDPIterator, typename TDPValue, typename TSeqVIterator,
//    typename TSeqHValue, typename TDPFormula, typename TDPPhase, typename TAlignProfile, typename TBandSpec>
//inline bool _computeFirstCell(TTraceTracker & tracker,
//                              TTraceIterator & traceIter,
//                              TDPIterator & dpIter,
//                              TDPValue const & prevD,
//                              TDPValue const & prevH,
//                              TDPValue const & prevV,
//                              TSeqVIterator const & seqVIter,
//                              TSeqVIterator const & /*seqVStop*/,
//                              TSeqHValue const & seqHValue,
//                              TDPFormula const & dpFormula,
//                              TDPPhase const & /*dpPhase*/,
//                              DPManager<TAlignProfile, Band<TBandSpec> > const & /*dpManager*/)
//{
//    typedef typename SetUpColumnManager<TAlignProfile, TDPPhase>::Type TColumnManager;
//    typedef Band<TBandSpec> TBand;
//    typedef typename IsTracebackOn<TAlignProfile>::Type TIsTracebackOn;
//
//    computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
//                 typename TColumnManager::TFirstCell(), TIsTracebackOn());
//    return false;
//}
//
//
//// ----------------------------------------------------------------------------
//// Function _computeLastCell
//// ----------------------------------------------------------------------------
//
//template <typename TTraceTracker, typename TTraceIterator, typename TDPIterator, typename TDPValue, typename TSeqVValue,
//    typename TSeqHValue, typename TDPFormula, typename TAlignProfile, typename TBandSpec>
//inline void _computeLastCell(TTraceTracker & tracker,
//                              TTraceIterator & traceIter,
//                              TDPIterator & dpIter,
//                              TDPValue const & prevD,
//                              TDPValue const & prevH,
//                              TDPValue const & prevV,
//                              TSeqVValue const & seqVValue,
//                              TSeqHValue const & seqHValue,
//                              bool const & isLastCell,
//                              TDPFormula const & dpFormula,
//                              DPBandClosurePhase const & /*dpPhase*/,
//                              DPManager<TAlignProfile, Band<TBandSpec> > const & /*dpManager*/)
//{
//    typedef typename SetUpColumnManager<TAlignProfile, DPBandClosurePhase>::Type TColumnManager;
//    typedef Band<TBandSpec> TBand;
//    typedef typename IsTracebackOn<TAlignProfile>::Type TIsTracebackOn;
//
//    if (!isLastCell)
//    {
//        computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, seqVValue, dpFormula,
//                     typename TColumnManager::TLastCell(), TIsTracebackOn());
//    }
//}
//
//template <typename TTraceTracker, typename TTraceIterator, typename TDPIterator, typename TDPValue, typename TSeqVValue,
//    typename TSeqHValue, typename TDPFormula, typename TAlignProfile, typename TBandSpec>
//inline void _computeLastCell(TTraceTracker & tracker,
//                              TTraceIterator & traceIter,
//                              TDPIterator & dpIter,
//                              TDPValue const & prevD,
//                              TDPValue const & prevH,
//                              TDPValue const & prevV,
//                              TSeqVValue const & seqVValue,
//                              TSeqHValue const & seqHValue,
//                              bool const & isLastCell,
//                              TDPFormula const & dpFormula,
//                              DPWideBandClosurePhase const & /*dpPhase*/,
//                              DPManager<TAlignProfile, Band<TBandSpec> > const & /*dpManager*/)
//{
//    typedef typename SetUpColumnManager<TAlignProfile, DPWideBandClosurePhase>::Type TColumnManager;
//    typedef Band<TBandSpec> TBand;
//    typedef typename IsTracebackOn<TAlignProfile>::Type TIsTracebackOn;
//
//    if (!isLastCell)
//    {
//        computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, seqVValue, dpFormula,
//                     typename TColumnManager::TLastCell(), TIsTracebackOn());
//    }
//}
//
//template <typename TTraceTracker, typename TTraceIterator, typename TDPIterator, typename TDPValue, typename TSeqVValue,
//    typename TSeqHValue, typename TDPFormula, typename TDPPhase, typename TAlignProfile, typename TBandSpec>
//inline void _computeLastCell(TTraceTracker & tracker,
//                              TTraceIterator & traceIter,
//                              TDPIterator & dpIter,
//                              TDPValue const & prevD,
//                              TDPValue const & prevH,
//                              TDPValue const & prevV,
//                              TSeqVValue const & seqVValue,
//                              TSeqHValue const & seqHValue,
//                              bool const & /*isLastCell*/,
//                              TDPFormula const & dpFormula,
//                              TDPPhase const & /*dpPhase*/,
//                              DPManager<TAlignProfile, Band<TBandSpec> > const & /*dpManager*/)
//{
//    typedef typename SetUpColumnManager<TAlignProfile, TDPPhase>::Type TColumnManager;
//    typedef Band<TBandSpec> TBand;
//    typedef typename IsTracebackOn<TAlignProfile>::Type TIsTracebackOn;
//
//    computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, seqVValue, dpFormula,
//                 typename TColumnManager::TLastCell(), TIsTracebackOn());
//}

// ----------------------------------------------------------------------------
// Function fillColumn()
// ----------------------------------------------------------------------------

template <typename TTraceTracker, typename TTraceIterator, typename TDPIterator, typename TSeqVIterator,
typename TSeqHValue, typename TDPFormula, typename TAlignProfile,  typename TBandSpec>
inline void fillColumn(TTraceTracker & tracker,
                       TTraceIterator & traceIter,
                       TDPIterator & dpIter,
                       TSeqVIterator & seqVIter,
                       TSeqVIterator const & seqVStop,
                       TSeqHValue const & seqHValue,
                       TDPFormula const & dpFormula,
                       DPWideBandMiddlePhase const & /*columnManger*/,
                       DPManager<TAlignProfile, Band<BandSwitchedOn<TBandSpec> > > const & /*dpManager*/)
{
        // computing the first cell in banded alignment
    typedef typename SetUpColumnManager<TAlignProfile, DPWideBandMiddlePhase>::Type TColumnManager;
    typedef typename IsTracebackOn<TAlignProfile>::Type TIsTracebackOn;
    typedef typename Value<TDPIterator>::Type TDPValue;
    typedef Band<BandSwitchedOff> TBand;

    TDPValue prevH = value(dpIter);
    TDPValue prevD;
    TDPValue prevV;
    _setInitValue(dpIter, tracker, TBand(), TColumnManager());
    computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
                 typename TColumnManager::TFirstCell(), TIsTracebackOn());
//    std::cout << value(dpIter) << "\t";
    ++seqVIter;
    increment(traceIter, TIsTracebackOn());
    while (seqVIter != seqVStop)
    {
        prevD = prevH;
        prevV = value(dpIter);
        prevH = value(++dpIter);
        computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
                     typename TColumnManager::TInnerCell(), TIsTracebackOn());
//        std::cout << value(dpIter) << "\t";
        ++seqVIter;
        increment(traceIter, TIsTracebackOn());
    }
    prevD = prevH;
    prevV = value(dpIter);
    prevH = value(++dpIter);
    computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
                 typename TColumnManager::TLastCell(), TIsTracebackOn());
//    std::cout << value(dpIter) << "\n";
}

template <typename TTraceTracker, typename TTraceIterator, typename TDPIterator, typename TSeqVIterator,
    typename TSeqHValue, typename TDPFormula, typename TDPPhase, typename TAlignProfile, typename TBandSpec>
inline void fillColumn(TTraceTracker & tracker,
                       TTraceIterator & traceIter,
                       TDPIterator & dpIter,
                       TSeqVIterator & seqVIter,
                       TSeqVIterator const & seqVStop,
                       TSeqHValue const & seqHValue,
                       TDPFormula const & dpFormula,
                       TDPPhase const & dpPhase,
                       DPManager<TAlignProfile, Band<BandSwitchedOn<TBandSpec> > > const & /*dpManager*/)
{
    // computing the first cell in banded alignment
    typedef typename SetUpColumnManager<TAlignProfile, TDPPhase>::Type TColumnManager;
    typedef Band<TBandSpec> TBand;
    typedef typename IsTracebackOn<TAlignProfile>::Type TIsTracebackOn;
    typedef typename Value<TDPIterator>::Type TDPValue;

    std::cout << "######## A" << std::endl;
    TDPIterator prevIter = dpIter;
    std::cout << "######## B" << std::endl;
    TDPValue prevH = value(++prevIter);
    std::cout << "######## C" << std::endl;
    TDPValue prevD = value(dpIter); // invalid read here since it does not point correctly ...
    std::cout << "######## D" << std::endl;
    TDPValue prevV;
    std::cout << "######## E" << std::endl;
    _setInitValue(dpIter, tracker, TBand(), dpPhase);  // only used within banded chain alignment to set the initialization value determined in the previous DP algorithm.
    std::cout << "######## F" << std::endl;
    computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
                 typename TColumnManager::TFirstCell(), TIsTracebackOn());
    std::cout << "######## G" << std::endl;
//    std::cout << value(dpIter) << "\t";
    if (seqVIter == seqVStop)  // TODO(rmaerker): What influence does this check has on the runtime?
    {
        return;
    }
    ++seqVIter;
//    std::cout << value(seqVIter) << " ";
    increment(traceIter, TIsTracebackOn());
    while (seqVIter != seqVStop)
    {
        prevD = prevH;
        prevV = value(dpIter);
        ++dpIter;
        prevH = value(++prevIter);
        computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
                     typename TColumnManager::TInnerCell(), TIsTracebackOn());
//        std::cout << value(dpIter) << "\t";
        increment(traceIter, TIsTracebackOn());
        ++seqVIter;
//        std::cout << value(seqVIter) << " ";
    }

    prevD = prevH;
    prevV = value(dpIter);
    ++dpIter;
    _getHorizontalValue(prevH, prevIter, dpPhase);
    computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
                 typename TColumnManager::TLastCell(), TIsTracebackOn());

//    std::cout << value(dpIter);
    // Note: that the dpIter always is == end(dpColumn)
}

template <typename TTraceTracker, typename TTraceIterator, typename TDPIterator, typename TSeqVIterator,
typename TSeqHValue, typename TDPFormula, typename TColumnType, typename TAlignProfile>
inline void fillColumn(TTraceTracker & tracker,
                       TTraceIterator & traceIter,
                       TDPIterator & dpIter,
                       TSeqVIterator & seqVIter,
                       TSeqVIterator const & seqVStop,
                       TSeqHValue const & seqHValue,
                       TDPFormula const & dpFormula,
                       TColumnType const & /*columnManger*/,
                       DPManager<TAlignProfile, Band<BandSwitchedOff> > const & /*dpManager*/)
{
        // computing the first cell in banded alignment
    typedef typename SetUpColumnManager<TAlignProfile, TColumnType>::Type TColumnManager;
    typedef typename IsTracebackOn<TAlignProfile>::Type TIsTracebackOn;
    typedef typename Value<TDPIterator>::Type TDPValue;
    typedef Band<BandSwitchedOff> TBand;

    TDPValue prevH = value(dpIter);
    TDPValue prevD;
    TDPValue prevV;
    _setInitValue(dpIter, tracker, TBand(), TColumnManager());
    computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
                 typename TColumnManager::TFirstCell(), TIsTracebackOn());
//    std::cout << value(dpIter) << "\t";
    increment(traceIter, TIsTracebackOn());
    while (seqVIter != seqVStop)
    {
        prevD = prevH;
        prevV = value(dpIter);
        prevH = value(++dpIter);
        computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
                     typename TColumnManager::TInnerCell(), TIsTracebackOn());
//        std::cout << value(dpIter) << "\t";
        ++seqVIter;
        increment(traceIter, TIsTracebackOn());
    }
    prevD = prevH;
    prevV = value(dpIter);
    prevH = value(++dpIter);
    computeCell(tracker, traceIter, value(dpIter), prevD, prevH, prevV, seqHValue, value(seqVIter), dpFormula,
                 typename TColumnManager::TLastCell(), TIsTracebackOn());
//    std::cout << value(dpIter) << "\n";
}

// ----------------------------------------------------------------------------
// function initializeMatrix
// ----------------------------------------------------------------------------

template<typename TTraceTracker, typename TTraceIterator, typename TDPIterator,
        typename TSeqVIterator, typename TSeqHIterator, typename TDPFormula, typename TDPPhase,
        typename TAlignmentProfile, typename TBand>
inline void initializeMatrix(TTraceTracker & tracker,
        TTraceIterator & traceIter,
        TDPIterator & activeColIter,
        TSeqVIterator const & seqVBegin,
        TSeqHIterator const & seqHIter,
        TDPFormula const & dpFormula,
        TDPPhase const & dpPhase,
        DPManager<TAlignmentProfile, TBand> & dpManager)
{
    TSeqVIterator seqVIter = seqVBegin + dpManager._spanSeqVBegin; // set begin of vertical seq relative to band
    TSeqVIterator seqVStop = seqVBegin + dpManager._spanSeqVEnd;  // set end of vertical sequence relative to band
    fillColumn(tracker, traceIter, activeColIter, seqVIter, seqVStop, value(seqHIter), dpFormula,
               dpPhase, dpManager);
    std::cout << std::endl;
}

// ----------------------------------------------------------------------------
// function fillMatrix
// ----------------------------------------------------------------------------

template<typename TTraceTracker, typename TTraceIterator, typename TDPIterator,
        typename TSeqVIterator, typename TSeqHIterator, typename TDPFormula, typename TDPPhase,
        typename TAlignmentProfile, typename TBand>
inline void fillMatrix(TTraceTracker & tracker,
        TTraceIterator & traceIter,
        TDPIterator & activeColIter,
        TSeqVIterator const & seqVBegin,
        TSeqHIterator & seqHIter,
        TSeqHIterator const & seqHStop,
        TDPFormula const & dpFormula,
        TDPPhase const & dpPhase,
        DPManager<TAlignmentProfile, TBand> & dpManager)
{
    // here set the value
    TSeqVIterator seqVIter;
    TSeqVIterator seqVStop;
    //go over seqHori

    for (; seqHIter != seqHStop; ++seqHIter)
    {
//        std::cout << value(seqHIter) << " - ";
        update(dpManager, dpPhase);
        seqVIter = seqVBegin + dpManager._spanSeqVBegin; // set begin of vertical seq relative to band
        seqVStop = seqVBegin + dpManager._spanSeqVEnd;  // set end of vertical sequence relative to band
        traceIter += dpManager._spanTrace;
        activeColIter += dpManager._spanDp;
        goBeginNextColumn(tracker, dpManager._spanSeqVBegin, TBand());  // used if tracker needs the correct positions in the grid.
        fillColumn(tracker, traceIter, activeColIter, seqVIter, seqVStop, value(seqHIter), dpFormula, dpPhase, dpManager);
//        std::cout << value(seqHIter) << "\n";
//        for (TDPIterator testIter = activeColIter - 13; testIter != activeColIter; ++testIter)
//        {
//            std::cout << value(testIter) << " ";
//        }
//        std::cout << std::endl;
    }
}

// ----------------------------------------------------------------------------
// function computeDPMatrix
// ----------------------------------------------------------------------------

template<typename TTraceTracker, typename TTraceMatrix, typename TSequenceH, typename TSequenceV, typename TScoreScheme,
        typename TBand, typename TAlignmentProfile>
void computeDPMatrix(TTraceTracker & tracker,
        TTraceMatrix & traceMatrix,
        TSequenceH const & seqH,
        TSequenceV const & seqV,
        TScoreScheme const & scoreScheme,
        TBand const & band,
        TAlignmentProfile const & /*alignmentProfile*/)
{
    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignSpec;
    typedef typename Value<TScoreScheme>::Type TScoreValue;
    typedef typename GetGapSpec<TAlignmentProfile>::Type TGapSpec;
    typedef typename Value<TScoreScheme>::Type TScoreValue;
    typedef DPValue<TScoreValue, TGapSpec> TDPValue;
    typedef DPMatrix<TDPValue> TDPMatrix;
    typedef typename Size<TDPMatrix>::Type TSize;
    typedef typename Position<TDPMatrix>::Type TPosition;
    typedef typename Iterator<TDPMatrix>::Type TDPIterator;

    typedef typename TraceBitMask::Type TTraceValue;
    typedef typename Iterator<TTraceMatrix>::Type TTraceIterator;
    typedef typename IsTracebackOn<TAlignmentProfile>::Type TTracebackOn;

    typedef typename Iterator<TSequenceH const>::Type TSeqHIterator;
    typedef typename Iterator<TSequenceV const>::Type TSeqVIterator;

    typedef DPManager<TAlignmentProfile, TBand> TDPManager;


    // Check valid DP settings

    if(!_isValidDPSettings(seqH, seqV, band, TAlignmentProfile()))
    {
        resize(tracker._tracePoints, 1);
        value(begin(traceMatrix)) = TraceBitMask::NONE;
        tracker._maxScore = MinValue<TScoreValue>::VALUE;
        tracker._tracePoints[0] = begin(traceMatrix);
        return;
    }

    // Parameter Setup

    // column parameter
    TSize columnSize = getColumnSize(seqH, seqV, band);
    TSize initialColumnSize = getInitialColumnSize(seqH, seqV, band);

    TSeqVIterator seqVBegin = begin(seqV);

    // DP Matrix
    TDPMatrix dpMatrix(seqH, seqV, band);  // define the dp matrix given both sequences

    // Column Manager
    TDPManager dpManager(columnSize, initialColumnSize);

    // DP Iterator
    TDPIterator dpIter = begin(dpMatrix, dpManager);

    // Traceback iterator
    TTraceIterator traceIter = begin(traceMatrix, dpManager);

    // DP Formula
    DPFormula<TScoreScheme, TAlignSpec> dpFormula(scoreScheme);

    std::cout << seqH << std::endl;
    // here we run the recursion but based on the band we set up the columns differently
    if (!_isBandEnabled(band))
    {
        // case: F, L, O
        //INITIALIZATION
        initializeMatrix(tracker, traceIter, dpIter, seqVBegin, begin(seqH), dpFormula, DPNoBandInitPhase(),
                         dpManager);

        //RECURSION
        TSeqHIterator seqHBegin = begin(seqH); // is begin of seqH
        TSeqHIterator seqHEnd = end(seqH);

        // second phase full band or no band
        fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHBegin, seqHEnd, dpFormula,
                DPNoBandPhase(), dpManager);
    }
    else if (!_isWideBand(band, seqH, seqV))
    {
        // case: A, C, D, E, G, I, J
        //INITIALIZATION
        // TODO(rmaerker): Check if we need DPNoBandInit...

        initializeMatrix(tracker, traceIter, dpIter, seqVBegin, begin(seqH), dpFormula, DPBandInitPhase(),
                         dpManager);

        //RECURSION

        TSeqHIterator seqHBeginFirstPhase = _beginBandFirstPhase(seqH, seqV, band); // is begin of seqH
        TSeqHIterator seqHBandedMiddlePhase = _beginBandMiddlePhase(seqH, seqV, band);
        TSeqHIterator seqHBandedLastPhase = _beginBandEndPhase(seqH, seqV, band);
        TSeqHIterator seqHEnd = _endBand(seqH, seqV, band);

        fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHBeginFirstPhase, seqHBandedMiddlePhase, dpFormula,
                DPBandFirstPhase(), dpManager);

        fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHBandedMiddlePhase, seqHBandedLastPhase, dpFormula,
                DPBandMiddlePhase(), dpManager);

        // not necessarily need to track this. For instance if we are not at the end of the sequence.
        track(tracker, value(dpIter), traceIter, typename IsFreeEndGap<TAlignmentProfile, LastRow>::Type());

        _correctSpan(dpManager);
        fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHBandedLastPhase, seqHEnd, dpFormula,
                DPBandLastPhase(), dpManager);
    }
    else
    { // case B, H, K
        //INITIALIZATION
        initializeMatrix(tracker, traceIter, dpIter, seqVBegin, begin(seqH), dpFormula, DPWideBandInitPhase(),
                         dpManager);
        //RECURSION

        TSeqHIterator seqHBeginFirstPhase = _beginWideBandFirstPhase(seqH, seqV, band); // is begin of seqH
        TSeqHIterator seqHBandedMiddlePhase = _beginWideBandMiddlePhase(seqH, seqV, band);
        TSeqHIterator seqHBandedLastPhase = _beginWideBandEndPhase(seqH, seqV, band);
        TSeqHIterator seqHEnd = _endWideBand(seqH, seqV, band);

        fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHBeginFirstPhase, seqHBandedMiddlePhase, dpFormula,
                DPWideBandFirstPhase(), dpManager);
//        std::cout << "After First " << seqHBandedMiddlePhase - begin(seqH) << std::endl;

        // TODO(rmaerker): how to track last cell if not enabled...
        track(tracker, value(dpIter), traceIter, typename IsFreeEndGap<TAlignmentProfile, LastRow>::Type());
//        std::cout << "After Intermed " << seqHBandedMiddlePhase - begin(seqH) << seqHBandedLastPhase1 - begin(seqH) << std::endl;
        fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHBandedMiddlePhase, seqHBandedLastPhase, dpFormula,
                DPWideBandMiddlePhase(), dpManager);

//        std::cout << "After Middle " << seqHBandedLastPhase - begin(seqH) << std::endl;

//        fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHBandedLastPhase, seqHBandedLastPhase2, dpFormula,
//                DPWideBandMiddlePhase(), dpManager);
//        std::cout << "After Middle " <<  seqHBandedLastPhase2 - begin(seqH) << std::endl;
        _correctSpan(dpManager);
        fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHBandedLastPhase, seqHEnd, dpFormula,
                DPWideBandLastPhase(), dpManager);
//        std::cout << "After Last " <<  seqHBandedLastPhase - begin(seqH) << " - " << seqHEnd - begin(seqH) << std::endl;
    }
    // TODO(rmaerker): special case of only one diagonal ....

    // We need to scan the last value
    trackLastColumn(tracker, dpIter, (TDPIterator) begin(dpMatrix), traceIter, band,
                    typename IsFreeEndGap<TAlignmentProfile, LastColumn>::Type());
//
//    //first phase initialization of the band
//    fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHIter, seqHStop1, dpFormula,
//            BandOpenColumn(), dpManager);
//    // second phase full band or no band
//    fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHIter, seqHStop2, dpFormula,
//            FullColumn(), dpManager);
//    // third phase end of band
//    fillMatrix(tracker, traceIter, dpIter, seqVBegin, seqHIter, seqHStop3, dpFormula,
//            BandCloseColumn(), dpManager);
    // trace back last column and track scores ... either last cell only or last column depending on alignment profile

}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_IMPL_H_
