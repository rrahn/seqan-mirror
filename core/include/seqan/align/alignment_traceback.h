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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_H_

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

//---------------------------------------------------
// function getTraceMatrixSize()
//---------------------------------------------------

template <typename TSequence0, typename TSequence1, typename TBand>
inline typename Size<TSequence0>::Type
getTraceMatrixSize(TSequence0 const & seq0,
                   TSequence1 const & seq1,
                   TBand const & band,
                   True const & /*traceback on*/)
{
    return getColumnSize(seq0, seq1, band) * (length(seq0)+1);
}

template <typename TSequence0, typename TSequence1, typename TBand>
inline typename Size<TSequence0>::Type
getTraceMatrixSize(TSequence0 const & /*seq0*/,
                   TSequence1 const & /*seq1*/,
                   TBand const & /*band*/,
                   False const & /*traceback off*/)
{   // return 1 so that begin of trace matrix is valid
    return typename Size<TSequence0>::Type(1);
}

// ----------------------------------------------------------------------------
// function computeTraceback()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TTracker, typename TTraceMatrix, typename TSequence0, typename TSequence1,
        typename TBand, typename TCore, typename TGaps>
inline void
computeTraceback(TTarget const & /*target*/,
                 TTracker const & /*tracker*/,
                 TTraceMatrix const & /*traceMatrix*/,
                 TSequence0 const & /*seq0*/,
                 TSequence1 const & /*seq1*/,
                 TBand const & /*band*/,
                 AlignmentProfile<TCore, TGaps, TracebackSwitchedOff> const & /*alignmentProfile*/)
{
    // nothing to do
}

template <typename TTraceSegment, typename TTracker, typename TTraceMatrix, typename TSequenceHorizontal,
        typename TSequenceVertical, typename TBand, typename TCore, typename TGaps>
inline void
computeTraceback(String<TTraceSegment> & target,
                 TTracker & tracker,
                 TTraceMatrix const & traceMatrix,
                 TSequenceHorizontal const & seqHorizontal,
                 TSequenceVertical const & seqVertical,
                 TBand const & band,
                 AlignmentProfile<TCore, TGaps, TracebackSwitchedOn> const & alignmentProfile)
{
    typedef typename TTracker::TTracePoints TTracePoints;
    typedef typename Iterator<TTracePoints>::Type TIterator;

//    printTraceback(traceMatrix, seqHorizontal, seqVertical, band);

    followTrace(target, tracker._tracePoints[0], traceMatrix, seqHorizontal, seqVertical, band, alignmentProfile);

}

template <typename TTraceSegment, typename TTracker, typename TTraceMatrix, typename TSequenceHorizontal,
        typename TSequenceVertical, typename TBand, typename TCore, typename TGaps>
inline void
computeTraceback(StringSet<String<TTraceSegment> > & target,
                 TTracker & tracker,
                 TTraceMatrix const & traceMatrix,
                 TSequenceHorizontal const & seqHorizontal,
                 TSequenceVertical const & seqVertical,
                 TBand const & band,
                 AlignmentProfile<TCore, TGaps, TracebackSwitchedOn> const & alignmentProfile)
{
    typedef typename TTracker::TTracePoints TTracePoints;
    typedef typename Iterator<TTracePoints>::Type TIterator;

    TIterator it = begin(tracker._tracePoints);
    TIterator itEnd = end(tracker._tracePoints);

    unsigned idx = 0;
    for (; it != itEnd; ++it, ++idx)
    {
        followTrace(target[idx], *it, traceMatrix, seqHorizontal, seqVertical, band, alignmentProfile);
    }

}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_H_
