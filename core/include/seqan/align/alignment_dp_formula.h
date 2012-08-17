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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_FORMULA_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_FORMULA_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================


template <typename TScoreScheme, typename TAlignmentSpec>
struct DPFormula{};

template <typename TScoreScheme, typename TGlobalSpec>
struct DPFormula<TScoreScheme, Global<TGlobalSpec> >
{
    typedef typename TraceBitMask::Type TTraceValue;

    TScoreScheme _scoreScheme;

    DPFormula() : _scoreScheme()
    {}

    DPFormula(TScoreScheme const & other) : _scoreScheme(other)
    {}

    template <typename TDPValue, typename TSeq0Value, typename TSeq1Value, typename TDPDirection>
    inline TTraceValue
    operator()(TDPValue & activeCell,
               TDPValue const & prevDiagonal,
               TDPValue const & prevHorizontal,
               TDPValue const & prevVertical,
               TSeq0Value const & seq0Value,
               TSeq1Value const & seq1Value,
               TDPDirection const & dpDirection) const
    {
        return computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, seq0Value, seq1Value,
                                         _scoreScheme, dpDirection);
    }

    template <typename TScoreValue, typename TDPValueSpec, typename TSeq0Value, typename TSeq1Value, typename TDPDirection>
    inline TTraceValue
    operator()(DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> & activeCell,
               DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> const & prevDiagonal,
               DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> const & prevHorizontal,
               DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> const & prevVertical,
               TSeq0Value const & seq0Value,
               TSeq1Value const & seq1Value,
               TDPDirection const & dpDirection) const
    {
        TTraceValue trace = computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, seq0Value, seq1Value,
                                         _scoreScheme, dpDirection);
        // additionally check if cell is forbidden
        if ((trace | TraceBitMask::FORBIDDEN) == TraceBitMask::FORBIDDEN)
        {
            setForbidden(activeCell, true);
        }
        return trace;
    }
};

template <typename TScoreScheme, typename TLocalSpec>
struct DPFormula<TScoreScheme, Local<TLocalSpec> >
{
    typedef typename TraceBitMask::Type TTraceValue;
    typedef typename Value<TScoreScheme>::Type TScoreValue;

    TScoreScheme _scoreScheme;

    DPFormula() : _scoreScheme(){ }

    DPFormula(TScoreScheme const & other) : _scoreScheme(other){ }



    template <typename TDPValue, typename TSeq0Value, typename TSeq1Value, typename TDPDirection>
    inline TTraceValue
    operator()(TDPValue & activeCell,
               TDPValue const & prevDiagonal,
               TDPValue const & prevHorizontal,
               TDPValue const & prevVertical,
               TSeq0Value const & seq0Value,
               TSeq1Value const & seq1Value,
               TDPDirection const & dpDirection) const
    {
        TTraceValue trace = computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, seq0Value, seq1Value,
                                         _scoreScheme, dpDirection);
        if (getScore(activeCell) <= 0)
        {
            return _doComputeScoreZero(activeCell);
        }
        return trace;
    }

    template <typename TScoreValue, typename TDPValueSpec, typename TSeq0Value, typename TSeq1Value, typename TDPDirection>
    inline TTraceValue
    operator()(DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> & activeCell,
               DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> const &  prevDiagonal,
               DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> const &  prevHorizontal,
               DPValue<TScoreValue, TDPValueSpec, DPValueConfigForbidden> const &  prevVertical,
               TSeq0Value const & seq0Value,
               TSeq1Value const & seq1Value,
               TDPDirection const & dpDirection) const
    {
        TTraceValue trace = computeScore(activeCell, prevDiagonal, prevHorizontal, prevVertical, seq0Value, seq1Value,
                                         _scoreScheme, dpDirection);

        if ((trace | TraceBitMask::FORBIDDEN) == TraceBitMask::FORBIDDEN)
        {
            setForbidden(activeCell, true);
            return trace;
        }
        if (getScore(activeCell) <= 0)
        {
            return _doComputeScoreZero(activeCell);
        }
        return trace;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// function _doComputeScoreDiagonal()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TDPValueSpec, typename TDPValueConfig, typename TSeq0Value,
        typename TSeq1Value, typename TScoreScheme>
inline TScoreValue
_doComputeScoreDiagonal(DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> & activeCell,
                        DPValue<TScoreValue, TDPValueSpec, TDPValueConfig> const & prevCell,
                        TSeq0Value const & seq0Value,
                        TSeq1Value const & seq1Value,
                        TScoreScheme const & scoreScheme)
{
    if (isForbidden(prevCell)) // we can not come from a forbidden cell
    {
        activeCell = MinValue<TScoreValue>::VALUE; // caution with bit overflow
        return TraceBitMask::FORBIDDEN;
    }
    activeCell._score = prevCell._score + score(scoreScheme, seq0Value, seq1Value);
    return TraceBitMask::DIAGONAL;
}

// ----------------------------------------------------------------------------
// function _doComputeScoreVertical()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TDPValueConfig, typename TScoreScheme>
inline typename TraceBitMask::Type
_doComputeScoreVertical(DPValue<TScoreValue, LinearGaps, TDPValueConfig> & activeCell,
                        DPValue<TScoreValue, LinearGaps, TDPValueConfig> const & prevCell,
                        TScoreScheme const & scoreScheme)
{
    if (isForbidden(prevCell)) // we can not come from a forbidden cell
    {
        activeCell = MinValue<TScoreValue>::VALUE; // caution with bit overflow
        return TraceBitMask::FORBIDDEN;
    }
    activeCell._score = prevCell._score + scoreGap(scoreScheme);
    return TraceBitMask::VERTICAL;
}

template <typename TScoreValue, typename TDPValueConfig, typename TScoreScheme>
inline typename TraceBitMask::Type
_doComputeScoreVertical(DPValue<TScoreValue, AffineGaps, TDPValueConfig> & activeCell,
                          DPValue<TScoreValue, AffineGaps, TDPValueConfig> const & prevCell,
                          TScoreScheme const & scoreScheme)
{
    if (isForbidden(prevCell)) // we can not come from a forbidden cell
    {
        activeCell._scoreVertical = MinValue<TScoreValue>::VALUE; // caution with bit overflow
        return TraceBitMask::FORBIDDEN;
    }

    if (prevCell._dirForbidden & ForbiddenDirection::VERTICAL) // cannot extend gap
    {
        activeCell._scoreVertical = prevCell._score + scoreGapOpenVertical(scoreScheme);
        return TraceBitMask::VERTICAL_OPEN;
    }
    else
    {
        TScoreValue scoreOpen = prevCell._score + scoreGapOpenVertical(scoreScheme);
        TScoreValue scoreExtend = prevCell._scoreVertical + scoreGapExtendVertical(scoreScheme);
        if (scoreOpen > scoreExtend)
        {
            activeCell._scoreVertical = scoreOpen;
            return TraceBitMask::VERTICAL_OPEN;
        }
        activeCell._scoreVertical = scoreExtend;
    }
    return TraceBitMask::VERTICAL;
}

// ----------------------------------------------------------------------------
// function _doComputeScoreHorizontal()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TDPValueConfig, typename TScoreScheme>
inline typename TraceBitMask::Type
_doComputeScoreHorizontal(DPValue<TScoreValue, LinearGaps, TDPValueConfig> & activeCell,
                          DPValue<TScoreValue, LinearGaps, TDPValueConfig> const & prevCell,
                          TScoreScheme const & scoreScheme)
{
    if (isForbidden(prevCell)) // we can not come from a forbidden cell
    {
        activeCell = MinValue<TScoreValue>::VALUE; // caution with bit overflow
        return TraceBitMask::FORBIDDEN;
    }

    activeCell._score = prevCell._score + scoreGap(scoreScheme);
    return TraceBitMask::HORIZONTAL;
}

template <typename TScoreValue, typename TDPValueConfig, typename TScoreScheme>
inline typename TraceBitMask::Type
_doComputeScoreHorizontal(DPValue<TScoreValue, AffineGaps, TDPValueConfig> & activeCell,
                          DPValue<TScoreValue, AffineGaps, TDPValueConfig> const & prevCell,
                          TScoreScheme const & scoreScheme)
{
    if (isForbidden(prevCell))  // we can not come from a forbidden cell
    {
        activeCell._scoreHorizontal = MinValue<TScoreValue>::VALUE;  // caution with bit overflow
        return TraceBitMask::FORBIDDEN;
    }

    if (prevCell._dirForbidden & ForbiddenDirection::HORIZONTAL)  // gap extension is forbidden
    {
        activeCell._scoreHorizontal = prevCell._score + scoreGapOpenHorizontal(scoreScheme);
    }
    else
    {
        TScoreValue scoreOpen = prevCell._score + scoreGapOpenHorizontal(scoreScheme);
        TScoreValue scoreExtend = prevCell._scoreHorizontal + scoreGapExtendHorizontal(scoreScheme);
        if (scoreOpen > scoreExtend)
        {
            activeCell._scoreHorizontal = scoreOpen;
            return TraceBitMask::HORIZONTAL_OPEN;
        }
        activeCell._scoreHorizontal = scoreExtend;
    }
    return TraceBitMask::HORIZONTAL;
}

// ----------------------------------------------------------------------------
// function _doComputeScoreZero()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TConfig>
inline typename TraceBitMask::Type
_doComputeScoreZero(DPValue<TScoreValue, AffineGaps, TConfig> & activeCell)
{
    activeCell._score = static_cast<TScoreValue>(0);
    activeCell._scoreHorizontal = static_cast<TScoreValue>(0);
    activeCell._scoreVertical = static_cast<TScoreValue>(0);
    return TraceBitMask::NONE;
}

template <typename TScoreValue, typename TConfig>
inline typename TraceBitMask::Type
_doComputeScoreZero(DPValue<TScoreValue, LinearGaps, TConfig> & activeCell)
{
    activeCell._score = static_cast<TScoreValue>(0);
    return TraceBitMask::NONE;
}

// ----------------------------------------------------------------------------
// Function computeScore()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TDPValueConfig, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(DPValue<TScoreValue, AffineGaps, TDPValueConfig> & activeCell,
             DPValue<TScoreValue, AffineGaps, TDPValueConfig> const &  prevDiagonal,
             DPValue<TScoreValue, AffineGaps, TDPValueConfig> const &  prevHorizontal,
             DPValue<TScoreValue, AffineGaps, TDPValueConfig> const &  prevVertical,
             TSeq0Value const & seq0Value,
             TSeq1Value const & seq1Value,
             TScoreScheme const & scoreScheme,
             DPDirectionAll const & /*recursionSpec*/)
{
    typedef typename TraceBitMask::Type TTraceValue;

    TTraceValue trace = _doComputeScoreDiagonal(activeCell, prevDiagonal, seq0Value, seq1Value, scoreScheme);
    TTraceValue trace2 = _doComputeScoreVertical(activeCell, prevVertical, scoreScheme);

    if (activeCell._scoreVertical > activeCell._score)
    {
        activeCell._score = activeCell._scoreVertical;
        trace = trace2; // optimal score from vertical
    }
    else if (activeCell._scoreVertical == activeCell._score || (trace2 & TraceBitMask::VERTICAL_OPEN))
    {
        trace |= trace2; // optimal score from vertical and diagonal
    }

    trace2 = _doComputeScoreHorizontal(activeCell, prevHorizontal, scoreScheme);
    activeCell._dirForbidden = ForbiddenDirection::NONE;

    if (activeCell._scoreHorizontal > activeCell._score)
    {
        activeCell._score = activeCell._scoreHorizontal;
        return trace2; // optimal score from vertical
    }
    else if (activeCell._scoreHorizontal == activeCell._score || (trace2 & TraceBitMask::HORIZONTAL_OPEN))
    {
        return trace |= trace2; // optimal score coming from multiple directions
    }

    return trace;
}

template <typename TScoreValue, typename TDPValueConfig, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(DPValue<TScoreValue, LinearGaps, TDPValueConfig> & activeCell,
             DPValue<TScoreValue, LinearGaps, TDPValueConfig> const &  prevDiagonal,
             DPValue<TScoreValue, LinearGaps, TDPValueConfig> const &  prevHorizontal,
             DPValue<TScoreValue, LinearGaps, TDPValueConfig> const &  prevVertical,
             TSeq0Value const & seq0Value,
             TSeq1Value const & seq1Value,
             TScoreScheme const & scoreScheme,
             DPDirectionAll const & /*recursionSpec*/)
{
    typedef typename TraceBitMask::Type TTraceValue;
    typedef DPValue<TScoreValue, LinearGaps, TDPValueConfig> TDPValue;

    TTraceValue trace = _doComputeScoreVertical(activeCell, prevVertical, scoreScheme);
    TScoreValue tmpScore = getScore(activeCell);
    TTraceValue trace2 = _doComputeScoreHorizontal(activeCell, prevHorizontal, scoreScheme);
    if (getScore(activeCell) > tmpScore)
    {
        tmpScore = getScore(activeCell);
        trace = trace2; // optimal score from horizontal
    }
    else if (tmpScore == getScore(activeCell))
    {
        trace |= trace2; // optimal score from vertical and horizontal
    }
    trace2 = _doComputeScoreDiagonal(activeCell, prevDiagonal, seq0Value, seq1Value, scoreScheme);

    if (tmpScore > getScore(activeCell))
    {
        activeCell = tmpScore;
        return trace;
    }
    else if (tmpScore == getScore(activeCell))
    {
        return trace | trace2; // optimal score from multiple directions
    }
    return trace2;
}

template <typename TScoreValue, typename TDPValueConfig, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(DPValue<TScoreValue, AffineGaps, TDPValueConfig> & activeCell,
             DPValue<TScoreValue, AffineGaps, TDPValueConfig> const &  prevDiagonal,
             DPValue<TScoreValue, AffineGaps, TDPValueConfig> const &  /*prevHorizontal*/,
             DPValue<TScoreValue, AffineGaps, TDPValueConfig> const &  prevVertical,
            TSeq0Value const & seq0Value,
            TSeq1Value const & seq1Value,
            TScoreScheme const & scoreScheme,
            DPDirectionLowerBand const & /*recursionSpec*/)
{
    typedef typename TraceBitMask::Type TTraceValue;

    TTraceValue trace = _doComputeScoreDiagonal(activeCell, prevDiagonal, seq0Value, seq1Value, scoreScheme);
    TTraceValue trace2 = _doComputeScoreVertical(activeCell, prevVertical, scoreScheme);
    activeCell._dirForbidden = ForbiddenDirection::NONE;
    if (activeCell._scoreVertical > activeCell._score)
    {
        activeCell = activeCell._scoreVertical;
        return trace2; // optimal score from vertical
    }
    else if (activeCell._scoreVertical == getScore(activeCell) || (trace2 & TraceBitMask::VERTICAL_OPEN))
    {
        return (trace2 | trace); // optimal score from vertical and diagonal score
    }
    return trace; // optimal value from diagonal
}

template <typename TScoreValue, typename TDPValueConfig, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(DPValue<TScoreValue, LinearGaps, TDPValueConfig> & activeCell,
             DPValue<TScoreValue, LinearGaps, TDPValueConfig> const &  prevDiagonal,
             DPValue<TScoreValue, LinearGaps, TDPValueConfig> const &  /*prevHorizontal*/,
             DPValue<TScoreValue, LinearGaps, TDPValueConfig> const &  prevVertical,
            TSeq0Value const & seq0Value,
            TSeq1Value const & seq1Value,
            TScoreScheme const & scoreScheme,
            DPDirectionLowerBand const & /*recursionSpec*/)
{
    typedef DPValue<TScoreValue, LinearGaps, TDPValueConfig> TDPValue;
    typedef typename TraceBitMask::Type TTraceValue;

    TTraceValue trace2 = _doComputeScoreVertical(activeCell, prevVertical, scoreScheme);
    TScoreValue tmpScore = getScore(activeCell);
    TTraceValue trace = _doComputeScoreDiagonal(activeCell, prevDiagonal, seq0Value, seq1Value, scoreScheme);
    if(tmpScore > getScore(activeCell))
    {
        activeCell = tmpScore;
        return trace2;        // optimal score from vertical direction
    }
    else if (tmpScore == getScore(activeCell))
    {
        return (trace | trace2); // optimal score from both diagonal an vertical direction
    }
    return trace; // optimal score from diagonal
}


template <typename TScoreValue, typename TDPValueConfig, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(DPValue<TScoreValue, AffineGaps, TDPValueConfig> & activeCell,
             DPValue<TScoreValue, AffineGaps, TDPValueConfig> const &  prevDiagonal,
             DPValue<TScoreValue, AffineGaps, TDPValueConfig> const &  prevHorizontal,
             DPValue<TScoreValue, AffineGaps, TDPValueConfig> const &  /*prevVertical*/,
             TSeq0Value const & seq0Value,
             TSeq1Value const & seq1Value,
             TScoreScheme const & scoreScheme,
             DPDirectionUpperBand const & /*recursionSpec*/)
{
    typedef typename TraceBitMask::Type TTraceValue;
    TTraceValue trace = _doComputeScoreDiagonal(activeCell, prevDiagonal, seq0Value, seq1Value, scoreScheme);
    TTraceValue trace2 = _doComputeScoreHorizontal(activeCell, prevHorizontal, scoreScheme);
    activeCell._dirForbidden = ForbiddenDirection::NONE;
    if (activeCell._scoreHorizontal > activeCell._score)
    {
        activeCell = activeCell._scoreHorizontal;
        return trace2; // optimal score from vertical
    }
    else if (activeCell._scoreHorizontal == getScore(activeCell) || (trace2 & TraceBitMask::HORIZONTAL_OPEN))
    {
        return (trace | trace2); // optimal score from horizontal and diagonal score
    }
    return trace; // optimal value from diagonal
}

template <typename TScoreValue, typename TDPValueConfig, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(DPValue<TScoreValue, LinearGaps, TDPValueConfig> & activeCell,
             DPValue<TScoreValue, LinearGaps, TDPValueConfig> const &  prevDiagonal,
             DPValue<TScoreValue, LinearGaps, TDPValueConfig> const &  prevHorizontal,
             DPValue<TScoreValue, LinearGaps, TDPValueConfig> const &  /*prevVertical*/,
             TSeq0Value const & seq0Value,
             TSeq1Value const & seq1Value,
             TScoreScheme const & scoreScheme,
             DPDirectionUpperBand const & /*recursionSpec*/)
{
    typedef DPValue<TScoreValue, LinearGaps, TDPValueConfig> TDPValue;
    typedef typename TraceBitMask::Type TTraceValue;

    TTraceValue trace2 = _doComputeScoreHorizontal(activeCell, prevHorizontal, scoreScheme);
    TScoreValue tmpScore = getScore(activeCell);
    TTraceValue trace = _doComputeScoreDiagonal(activeCell, prevDiagonal, seq0Value, seq1Value, scoreScheme);

    if(tmpScore > getScore(activeCell))
    {
        activeCell = tmpScore;
        return trace2;    // optimal score from horizontal direction
    }
    else if (tmpScore == getScore(activeCell))
    {
        return (trace | trace2); // optimal score from diagonal and horizontal direction
    }
    return trace; // optimal value from diagonal
}

template <typename TDPValue, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(TDPValue & activeCell,
             TDPValue const &  prevDiagonal,
             TDPValue const &  /*prevHorizontal*/,
             TDPValue const &  /*prevVertical*/,
             TSeq0Value const & seq0Value,
             TSeq1Value const & seq1Value,
             TScoreScheme const & scoreScheme,
             DPDirectionDiagonal const & /*recursionSpec*/)
{
    setForbidden(activeCell, (ForbiddenDirection::VERTICAL | ForbiddenDirection::HORIZONTAL));
    return _doComputeScoreDiagonal(activeCell, prevDiagonal, seq0Value, seq1Value, scoreScheme);
}

template <typename TScoreValue, typename TConfig, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(DPValue<TScoreValue, AffineGaps, TConfig> & activeCell,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  prevDiagonal,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  /*prevHorizontal*/,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  /*prevVertical*/,
             TSeq0Value const & seq0Value,
             TSeq1Value const & seq1Value,
             TScoreScheme const & scoreScheme,
             DPDirectionDiagonal const & /*recursionSpec*/)
{
    activeCell._dirForbidden = (ForbiddenDirection::VERTICAL | ForbiddenDirection::HORIZONTAL);
    return _doComputeScoreDiagonal(activeCell, prevDiagonal, seq0Value, seq1Value, scoreScheme);
}

// used for vertical initialization
template <typename TDPValue, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(TDPValue & activeCell,
             TDPValue const &  /*prevDiagonal*/,
             TDPValue const &  /*prevHorizontal*/,
             TDPValue const &  prevVertical,
             TSeq0Value const & /*seq0Value*/,
             TSeq1Value const & /*seq1Value*/,
             TScoreScheme const & scoreScheme,
             DPDirectionVertical const & /*recursionSpec*/)
{
    return _doComputeScoreVertical(activeCell, prevVertical, scoreScheme);
}

template <typename TScoreValue, typename TConfig, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(DPValue<TScoreValue, AffineGaps, TConfig> & activeCell,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  /*prevDiagonal*/,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  /*prevHorizontal*/,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  prevVertical,
             TSeq0Value const & /*seq0Value*/,
             TSeq1Value const & /*seq1Value*/,
             TScoreScheme const & scoreScheme,
             DPDirectionVertical const & /*recursionSpec*/)
{
    typedef typename TraceBitMask::Type TTraceValue;
    TTraceValue trace = _doComputeScoreVertical(activeCell, prevVertical, scoreScheme);
    activeCell._dirForbidden = ForbiddenDirection::HORIZONTAL;
    activeCell._score = activeCell._scoreVertical;
    return trace;
}

// used for horizontal initialization

template <typename TDPValue, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(TDPValue & activeCell,
             TDPValue const &  /*prevDiagonal*/,
             TDPValue const &  prevHorizontal,
             TDPValue const &  /*prevVertical*/,
             TSeq0Value const & /*seq0Value*/,
             TSeq1Value const & /*seq1Value*/,
             TScoreScheme const & scoreScheme,
             DPDirectionHorizontal const & /*recursionSpec*/)
{
    return _doComputeScoreHorizontal(activeCell, prevHorizontal, scoreScheme);
}

template <typename TScoreValue, typename TConfig, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(DPValue<TScoreValue, AffineGaps, TConfig> & activeCell,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  /*prevDiagonal*/,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  prevHorizontal,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  /*prevVertical*/,
             TSeq0Value const & /*seq0Value*/,
             TSeq1Value const & /*seq1Value*/,
             TScoreScheme const & scoreScheme,
             DPDirectionHorizontal const & /*recursionSpec*/)
{
    typedef typename TraceBitMask::Type TTraceValue;
    TTraceValue trace = _doComputeScoreHorizontal(activeCell, prevHorizontal, scoreScheme);
    activeCell._dirForbidden = ForbiddenDirection::VERTICAL;
    activeCell._score = activeCell._scoreHorizontal;
    return trace;
}


template <typename TDPValue, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(TDPValue & activeCell,
             TDPValue const &  /*prevDiagonal*/,
             TDPValue const &  /*prevHorizontal*/,
             TDPValue const &  /*prevVertical*/,
             TSeq0Value const & /*seq0Value*/,
             TSeq1Value const & /*seq1Value*/,
             TScoreScheme const & /*scoreScheme*/,
             DPDirectionZero const & /*recursionSpec*/)
{
    return _doComputeScoreZero(activeCell);
}

template <typename TScoreValue, typename TConfig, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme>
inline typename TraceBitMask::Type
computeScore(DPValue<TScoreValue, AffineGaps, TConfig> & activeCell,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  /*prevDiagonal*/,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  /*prevHorizontal*/,
             DPValue<TScoreValue, AffineGaps, TConfig> const &  /*prevVertical*/,
             TSeq0Value const & /*seq0Value*/,
             TSeq1Value const & /*seq1Value*/,
             TScoreScheme const & /*scoreScheme*/,
             DPDirectionZero const & /*recursionSpec*/)
{
    activeCell._dirForbidden = (ForbiddenDirection::HORIZONTAL | ForbiddenDirection::VERTICAL);
    return _doComputeScoreZero(activeCell);
}

template <typename TDPValue, typename TSeq0Value, typename TSeq1Value, typename TScoreScheme, typename TBool>
inline typename TraceBitMask::Type
computeScore(TDPValue & activeCell,
             TDPValue const &  /*prevDiagonal*/,
             TDPValue const &  /*prevHorizontal*/,
             TDPValue const &  /*prevVertical*/,
             TSeq0Value const & /*seq0Value*/,
             TSeq1Value const & /*seq1Value*/,
             TScoreScheme const & /*scoreScheme*/,
             DPDirectionAdopt const & /*recursionSpec*/,
             TBool const & /*forbiddenActivated*/)
{
    if (isForbidden(activeCell))
         return TraceBitMask::FORBIDDEN;
    return TraceBitMask::NONE;
}

} // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_FORMULA_H_
