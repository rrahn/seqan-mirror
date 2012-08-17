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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_IMPL_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_IMPL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================


// ============================================================================
// Tags, Classes, Enums
// ============================================================================


//// TODO (rmaerker): remove debug code
//template <typename TTraceValue>
//String<char> translateTraceValue(TTraceValue const & traceValue, Simple const &)
//{
//    String<char> trace;
//    switch (traceValue)
//    {
//    case 0:
//        return "D";
//    case 1:
//        return "H";
//    case 2:
//        return "V";
//    case 3:
//        return "HV";
//    case 4:
//        return "Dh";
//    case 5:
//        return "Hh";
//    case 6:
//        return "Vh";
//    case 7:
//        return "HVh";
//    case 8:
//        return "Dv";
//    case 9:
//        return "Hv";
//    case 10:
//        return "Vv";
//    case 11:
//        return "HVv";
//    case 12:
//        return "Dhv";
//    case 13:
//        return "Hhv";
//    case 14:
//        return "Vhv";
//    case 15:
//        return "HVhv";
//    default:
//        return "X";
//    }
//
//    if (traceValue == 0)
//    {
//        appendValue(trace,"D");
//    }
//    else if (traceValue == 1)
//    {
//        return "H";
//    }
//    else if (traceValue == 2)
//    {
//        return "V";
//    }
//    return "X";
//}

//template <typename TMatrix, typename TSeqH, typename TSeqV, typename TBand>
//void printTraceback(TMatrix const & matrix, TSeqH const & seqH, TSeqV const & seqV, TBand const & band)
//{
//    typedef typename Iterator<TMatrix>::Type TIterator;
//    typedef typename Size<TSeqH>::Type TSize;
//
//    TSize columnSize = getColumnSize(seqH, seqV, band);
//    TSize rowSize = getRowSize(seqH, seqV, band);
//
//    std::cout << "\t\t";
//    for (TSize t = 0; t < length(seqH); ++t)
//    {
//        std::cout << seqH[t] << "\t";
//    }
//    std::cout << "\n";
//
//    for (TSize j = 0; j < columnSize; ++j)
//    {
//        if (j > 0)
//            std::cout << seqV[j-1] << "\t";
//        else
//            std::cout << "\t";
//        for (TSize i = 0; i < rowSize; ++i)
//        {
//            std::cout << translateTraceValue(matrix[i*columnSize + j]) << "\t";
//        }
//        std::cout << "\n";
//    }
//}

//template <typename TMatrix, typename TSeq>
//void printTraceback(TMatrix const & matrix, TSeq const & seqH, TSeq const & seqV)
//{
//    typedef typename Iterator<TMatrix>::Type TIterator;
//    typedef typename Size<TSeq>::Type TSize;
//
//    TSize columnSize = getColumnSize(seqH, seqV, Band<int, Off>());
//    TSize rowSize = getRowSize(seqH, seqV, Band<int, Off>());
//
//    std::cout << "\t\t";
//    for (TSize t = 0; t < length(seqH); ++t)
//    {
//        std::cout << seqH[t] << "\t";
//    }
//    std::cout << "\n";
//
//    for (TSize j = 0; j < columnSize; ++j)
//    {
//        if (j > 0)
//            std::cout << seqV[j-1] << "\t";
//        else
//            std::cout << "\t";
//        for (TSize i = 0; i < rowSize; ++i)
//        {
//            std::cout << translateTraceValue(matrix[i*columnSize + j], Simple()) << "\t";
//        }
//        std::cout << "\n";
//    }
//}

// ----------------------------------------------------------------------------
// function _goDiagonal()
// ----------------------------------------------------------------------------

template <typename TIterator, typename TSize>
inline void
_goDiagonal(TIterator & iter, TSize const & columnSize, Band<BandSwitchedOff> const & /*band*/)
{
    iter -= (columnSize + 1);
}

template <typename TIterator, typename TSize, typename TBandSpec>
inline void
_goDiagonal(TIterator & iter, TSize const & columnSize, Band<TBandSpec> const & /*band*/)
{
    iter -= columnSize;
}

// ----------------------------------------------------------------------------
// function _goVertial()
// ----------------------------------------------------------------------------

template <typename TIterator, typename TSize, typename TBandSpec>
inline void
_goVertical(TIterator & iter, TSize const & /*columnSize*/, Band<TBandSpec> const & /*band*/)
{
    --iter;
}

// ----------------------------------------------------------------------------
// function _goHorizontal()
// ----------------------------------------------------------------------------

template <typename TIterator, typename TSize>
inline void
_goHorziontal(TIterator & iter, TSize const & columnSize, Band<BandSwitchedOff> const & /*band*/)
{
    iter -= columnSize;
}

template <typename TIterator, typename TSize, typename TBandSpec>
inline void
_goHorziontal(TIterator & iter, TSize const & columnSize, Band<TBandSpec> const & /*band*/)
{
    iter -= (columnSize-1);
}

// ----------------------------------------------------------------------------
// function _determineSeq0Pos()
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSequenceH, typename TSequenceV, typename TBand>
inline TPosition _determineSeq0Pos(TPosition const & tracePosition,
                                   TSequenceH const & seqH,
                                   TSequenceV const & seqV,
                                   TBand const & band)
{
    return     tracePosition / getColumnSize(seqH, seqV, band);
}

// ----------------------------------------------------------------------------
// function _determineSeq1Pos()
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSequenceH, typename TSequenceV>
inline TPosition _determineSeq1Pos(TPosition const & tracePosition,
                                  TSequenceH const & seqH,
                                  TSequenceV const & seqV,
                                  Band<BandSwitchedOff> const & band)
{
    return tracePosition % getColumnSize(seqH, seqV, band);
}

template <typename TPosition, typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline TPosition _determineSeq1Pos(TPosition const & tracePosition,
                                  TSequenceH const & seqH,
                                  TSequenceV const & seqV,
                                  Band<TBandSpec> const & band)
{   // TODO (rmaerker): fix the function for upper levels
    typedef typename MakeSigned<TPosition>::Type TSignedPosition;
    TPosition currColumn = _determineSeq0Pos(tracePosition, seqH, seqV, band);
    // if currCol < getBandSize(band); then we need to remove the delta from the current position
    return _max(0,static_cast<TSignedPosition>(currColumn) - static_cast<TSignedPosition>(getUpperDiagonal(band)))
            + tracePosition % getColumnSize(seqH, seqV, band);
}

// ----------------------------------------------------------------------------
// Function _doFollowDiagonal()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TTraceIter, typename TTraceValue, typename TSeqHPos, typename TSeqVPos,
        typename TSize, typename TColumnSize, typename TBand, typename TGapSpec>
inline TTraceValue
_doFollowDiagonal(TTarget & target,
                  TTraceIter & traceIter,
                  TTraceValue & lastTraceValue,
                  TSeqHPos & seqHPos,
                  TSeqVPos & seqVPos,
                  TSize & fragmentLength,
                  TColumnSize const & columnSize,
                  TBand const & band,
                  TGapSpec const & /*tag*/)
{
    if (!(lastTraceValue & TraceBitMask::DIAGONAL)) // the old trace value was not diagonal
    {
        recordSegment(target, seqHPos, seqVPos, fragmentLength, lastTraceValue);

        lastTraceValue = TraceBitMask::DIAGONAL;
        fragmentLength = 0;
    }
    --seqHPos;
    --seqVPos;
    _goDiagonal(traceIter, columnSize, band);
    ++fragmentLength;
    return *traceIter;
}

// ----------------------------------------------------------------------------
// Function _doFollowVertical
// ----------------------------------------------------------------------------

template <typename TTarget, typename TTraceIter, typename TTraceValue, typename TSeqHPos, typename TSeqVPos,
        typename TSize, typename TColumnSize, typename TBand>
inline TTraceValue
_doFollowVertical(TTarget & target,
                  TTraceIter & traceIter,
                  TTraceValue & lastTraceValue,
                  TSeqHPos & seqHPos,
                  TSeqVPos & seqVPos,
                  TSize & fragmentLength,
                  TColumnSize const & columnSize,
                  TBand const & band,
                  AffineGaps const & /*tag*/)
{
    _doFollowVertical(target, traceIter, lastTraceValue, seqHPos, seqVPos, fragmentLength, columnSize, band,
                      LinearGaps());
    if (*traceIter & TraceBitMask::VERTICAL_OPEN)
    {
        --seqVPos;
        _goVertical(traceIter, columnSize, band);
        ++fragmentLength;
        return *traceIter;
    }
    return TraceBitMask::VERTICAL;
}

template <typename TTarget, typename TTraceIter, typename TTraceValue, typename TSeqHPos, typename TSeqVPos,
        typename TSize, typename TColumnSize, typename TBand>
inline TTraceValue
_doFollowVertical(TTarget & target,
                  TTraceIter & traceIter,
                  TTraceValue & lastTraceValue,
                  TSeqHPos & seqHPos,
                  TSeqVPos & seqVPos,
                  TSize & fragmentLength,
                  TColumnSize const & columnSize,
                  TBand const & band,
                  LinearGaps const & /*tag*/)
{
    if (!(lastTraceValue & TraceBitMask::VERTICAL))
    {
        recordSegment(target, seqHPos, seqVPos, fragmentLength, lastTraceValue);
        lastTraceValue = TraceBitMask::VERTICAL;
        fragmentLength = 0;
    }
    --seqVPos;
    _goVertical(traceIter, columnSize, band);
    ++fragmentLength;
    return *traceIter;
}

// ----------------------------------------------------------------------------
// Function _doFollowHorizontal()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TTraceIter, typename TTraceValue, typename TSeqHPos, typename TSeqVPos,
        typename TSize, typename TColumnSize, typename TBand>
inline TTraceValue
_doFollowHorizontal(TTarget & target,
                  TTraceIter & traceIter,
                  TTraceValue & lastTraceValue,
                  TSeqHPos & seqHPos,
                  TSeqVPos & seqVPos,
                  TSize & fragmentLength,
                  TColumnSize const & columnSize,
                  TBand const & band,
                  AffineGaps const & /*tag*/)
{
    _doFollowHorizontal(target, traceIter, lastTraceValue, seqHPos, seqVPos, fragmentLength, columnSize, band,
                      LinearGaps());
    if (*traceIter & TraceBitMask::HORIZONTAL_OPEN)
    { // follow one more
        --seqHPos;
        _goHorziontal(traceIter, columnSize, band);
        ++fragmentLength;
        return *traceIter;
    }
    return TraceBitMask::HORIZONTAL;
}

template <typename TTarget, typename TTraceIter, typename TTraceValue, typename TSeqHPos, typename TSeqVPos,
        typename TSize, typename TColumnSize, typename TBand>
inline TTraceValue
_doFollowHorizontal(TTarget & target,
                  TTraceIter & traceIter,
                  TTraceValue & lastTraceValue,
                  TSeqHPos & seqHPos,
                  TSeqVPos & seqVPos,
                  TSize & fragmentLength,
                  TColumnSize const & columnSize,
                  TBand const & band,
                  LinearGaps const & /*tag*/)
{
    if (!(lastTraceValue & TraceBitMask::HORIZONTAL))
    {
        recordSegment(target, seqHPos, seqVPos, fragmentLength, lastTraceValue);
        lastTraceValue = TraceBitMask::HORIZONTAL;
        fragmentLength = 0;
    }
    --seqHPos;
    _goHorziontal(traceIter, columnSize, band);
    ++fragmentLength;
    return *traceIter;
}

// ----------------------------------------------------------------------------
// Function _initializeTraceBack()
// ----------------------------------------------------------------------------

template <typename TTraceValue>
inline TTraceValue
_intializeTraceBack(TTraceValue & traceValue, AffineGaps const & /*tag*/)
{
    if (traceValue & TraceBitMask::VERTICAL)
    {
        traceValue &= (TraceBitMask::VERTICAL | TraceBitMask::VERTICAL_OPEN);
        return TraceBitMask::VERTICAL;
    }
    else if (traceValue & TraceBitMask::HORIZONTAL)
    {
        traceValue &= (TraceBitMask::HORIZONTAL | TraceBitMask::HORIZONTAL_OPEN);
        return TraceBitMask::HORIZONTAL;
    }
    else if (traceValue & TraceBitMask::DIAGONAL)
    {
        return TraceBitMask::DIAGONAL;
    }
    else if (traceValue & TraceBitMask::VERTICAL_OPEN)
    {
        traceValue = TraceBitMask::VERTICAL_OPEN;
        return TraceBitMask::VERTICAL;
    }
    else if (traceValue & TraceBitMask::HORIZONTAL_OPEN)
    {
        traceValue = TraceBitMask::HORIZONTAL_OPEN;
        return TraceBitMask::HORIZONTAL;
    }
    return TraceBitMask::NONE;
}

template <typename TTraceValue>
inline TTraceValue
_intializeTraceBack(TTraceValue const & traceValue, LinearGaps const & /*tag*/)
{
    if (traceValue & TraceBitMask::DIAGONAL)
    {
        return TraceBitMask::DIAGONAL;
    }
    else if (traceValue & TraceBitMask::VERTICAL)
    {
        return TraceBitMask::VERTICAL;
    }
    else if (traceValue & TraceBitMask::HORIZONTAL)
    {
        return TraceBitMask::HORIZONTAL;
    }
    return TraceBitMask::NONE;
}

// ----------------------------------------------------------------------------
// function _doFollowTrace()
// ----------------------------------------------------------------------------
template <typename TTarget, typename TTraceIterator, typename TSeq0Pos, typename TSeq1Pos, typename TSize,
    typename TBand, typename TGapSpec>
inline void
_doFollowTrace(TTarget & target,
               TTraceIterator & traceIter,
               TSeq0Pos & seqHPos,
               TSeq1Pos & seqVPos,
               TSize const & columnSize,
               TBand const & band,
               TGapSpec const & gapSpec)
{
    typedef typename TraceBitMask::Type TTraceValue;

    // the segement length value
    TSize fragmentLength = 0;
    TTraceValue traceValue = value(traceIter);
    TTraceValue lastTraceValue = _intializeTraceBack(traceValue, gapSpec);

    // TRACEBACK
    while ((seqHPos != 0 && seqVPos != 0) && traceValue != TraceBitMask::NONE)
    {

        if (traceValue & TraceBitMask::DIAGONAL)
        {
            traceValue = _doFollowDiagonal(target, traceIter, lastTraceValue, seqHPos, seqVPos, fragmentLength,
                                           columnSize, band, gapSpec);
        }
        else if ((traceValue & TraceBitMask::VERTICAL))
        {
            traceValue = _doFollowVertical(target, traceIter, lastTraceValue, seqHPos, seqVPos, fragmentLength,
                                           columnSize, band, gapSpec);
        }
        else if (traceValue & TraceBitMask::HORIZONTAL)
        {
            traceValue = _doFollowHorizontal(target, traceIter, lastTraceValue, seqHPos, seqVPos, fragmentLength,
                                           columnSize, band, gapSpec);
        }
        else if (traceValue & TraceBitMask::VERTICAL_OPEN)
        {
            traceValue = _doFollowVertical(target, traceIter, lastTraceValue, seqHPos, seqVPos, fragmentLength,
                                           columnSize, band, LinearGaps());
        }
        else if (traceValue & TraceBitMask::HORIZONTAL_OPEN)
        {
            traceValue = _doFollowHorizontal(target, traceIter, lastTraceValue, seqHPos, seqVPos, fragmentLength,
                                           columnSize, band, LinearGaps());
        }
        else
        { // the trace back is either NONE or something else
            if (traceValue == TraceBitMask::NONE)
            {
                break;
            }
            else if (traceValue & TraceBitMask::FORBIDDEN)
            {
                SEQAN_ASSERT_FAIL("Fatal Error in Traceback: Reached cell that was forbidden!");
            }
        }
    }
    // trace the last detected segment
    recordSegment(target, seqHPos, seqVPos, fragmentLength, lastTraceValue);
}

// ----------------------------------------------------------------------------
// function followTrace()
// ----------------------------------------------------------------------------

template<typename TTarget, typename TTraceIterator, typename TTraceMatrix, typename TSequenceH, typename TSequenceV,
        typename TBand, typename TCoreSpec, typename TGapSpec>
inline void followTrace(TTarget & target,
                        TTraceIterator & traceIter,
                        TTraceMatrix const & traceMatrix,
                        TSequenceH const & seqH,
                        TSequenceV const & seqV,
                        TBand const & band,
                        AlignmentProfile<Global<TCoreSpec>, TGapSpec, TracebackSwitchedOn> const &)
{
    typedef AlignmentProfile<Global<TCoreSpec>, TGapSpec, Traceback<Default> > TAlignmentProfile;
    typedef typename Size<TTraceMatrix>::Type TSize;
    typedef typename TraceBitMask::Type TTraceValue;
    typedef typename Position<TSequenceH>::Type TSeqHPos;
    typedef typename Position<TSequenceV>::Type TSeqVPos;


    // INITIALIZATION
    TTraceIterator traceItBegin = begin(traceMatrix);

    // Set the correct sequence positions
    TSeqHPos seqHPos = _determineSeq0Pos(traceIter - traceItBegin, seqH, seqV, band);
    TSeqVPos seqVPos = _determineSeq1Pos(traceIter - traceItBegin, seqH, seqV, band);

    // record overlaps at tail of sequences if any
    if (seqHPos != length(seqH))
    {
        recordSegment(target, seqHPos, seqVPos, length(seqH) - seqHPos, +TraceBitMask::HORIZONTAL);
    }
    if (seqVPos != length(seqV))
    {
        recordSegment(target, seqHPos, seqVPos, length(seqV) - seqVPos, +TraceBitMask::VERTICAL);
    }

    // actually follow the trace
    _doFollowTrace(target, traceIter, seqHPos, seqVPos, getColumnSize(seqH, seqV, band), band, TGapSpec());

    // record overlaps at begin of sequences if any
    if(seqHPos != 0)
    {
        recordSegment(target, 0, 0, seqHPos, +TraceBitMask::HORIZONTAL);
    }
    if(seqVPos != 0)
    {
        recordSegment(target, 0, 0, seqVPos, +TraceBitMask::VERTICAL);
    }
}


template<typename TTarget, typename TTraceIterator, typename TTraceMatrix, typename TSequenceH, typename TSequenceV,
        typename TBand, typename TCoreSpec, typename TGapSpec>
inline void followTrace(TTarget & target,
                        TTraceIterator & traceIter,
                        TTraceMatrix const & traceMatrix,
                        TSequenceH const & seqH,
                        TSequenceV const & seqV,
                        TBand const & band,
                        AlignmentProfile<Local<TCoreSpec>, TGapSpec, TracebackSwitchedOn> const &)
{
    typedef AlignmentProfile<Global<TCoreSpec>, TGapSpec, Traceback<Default> > TAlignmentProfile;
    typedef typename Size<TTraceMatrix>::Type TSize;
    typedef typename TraceBitMask::Type TTraceValue;
    typedef typename Position<TSequenceH>::Type TSeqHPos;
    typedef typename Position<TSequenceV>::Type TSeqVPos;


    // INITIALIZATION
    TTraceIterator traceItBegin = begin(traceMatrix);

    // Set the correct sequence positions
    TSeqHPos seqHPos = _determineSeq0Pos(traceIter - traceItBegin, seqH, seqV, band);
    TSeqVPos seqVPos = _determineSeq1Pos(traceIter - traceItBegin, seqH, seqV, band);

    // Note we don't record overlaps at the tail or at the beginning of the sequences.
    // Instead we cut the gaps by setting clip positions afterwards

    // actually follow the trace
    _doFollowTrace(target, traceIter, seqHPos, seqVPos, getColumnSize(seqH, seqV, band), band, TGapSpec());
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_IMPL_H_
