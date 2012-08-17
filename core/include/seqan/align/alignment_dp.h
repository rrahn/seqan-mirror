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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_H_

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

// ----------------------------------------------------------------------------
// SubstituteAlignConfig
// ----------------------------------------------------------------------------

template <typename TAlignConfig>
struct SubstituteAlignConfig;

// 0000
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<false, false, false, false, TSpec> >
{
    typedef FreeEndGaps<False, False, False, False> Type;
};

// 0001
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<false, false, false, true, TSpec> >
{
    typedef FreeEndGaps<False, False, True, False> Type;
};

// 0010
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<false, false, true, false, TSpec> >
{
    typedef FreeEndGaps<False, False, False, True> Type;
};


// 0011
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<false, false, true, true, TSpec> >
{
    typedef FreeEndGaps<False, False, True, True> Type;
};


// 0100
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<false, true, false, false, TSpec> >
{
    typedef FreeEndGaps<False, True, False, False> Type;
};


// 0101
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<false, true, false, true, TSpec> >
{
    typedef FreeEndGaps<False, True, True, False> Type;
};


// 0110
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<false, true, true, false, TSpec> >
{
    typedef FreeEndGaps<False, True, False, True> Type;
};


// 0111
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<false, true, true, true, TSpec> >
{
    typedef FreeEndGaps<False, True, True, True> Type;
};


// 1000
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<true, false, false, false, TSpec> >
{
    typedef FreeEndGaps<True, False, False, False> Type;
};


// 1001
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<true, false, false, true, TSpec> >
{
    typedef FreeEndGaps<True, False, True, False> Type;
};


// 1010
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<true, false, true, false, TSpec> >
{
    typedef FreeEndGaps<True, False, False, True> Type;
};


// 1011
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<true, false, true, true, TSpec> >
{
    typedef FreeEndGaps<True, False, True, True> Type;
};


// 1100
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<true, true, false, false, TSpec> >
{
    typedef FreeEndGaps<True, True, False, False> Type;
};


// 1101
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<true, true, false, true, TSpec> >
{
    typedef FreeEndGaps<True, True, True, False> Type;
};


// 1110
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<true, true, true, false, TSpec> >
{
    typedef FreeEndGaps<True, True, False, True> Type;
};

// 1111
template <typename TSpec>
struct SubstituteAlignConfig<AlignConfig<true, true, true, true, TSpec> >
{
    typedef FreeEndGaps<True, True, True, True> Type;
};

// ----------------------------------------------------------------------------
// SetUpAlignmentProfile
// ----------------------------------------------------------------------------

template <typename TAlgoTag, typename TAlignConfig, typename TTraceSwitch, typename TGapCosts>
struct SetUpAlignmentProfile;

// profile for NeedlemanWunsch
template <typename TAlignConfig, typename TGapCosts, typename TTraceSwitch>
struct SetUpAlignmentProfile<NeedlemanWunsch, TAlignConfig, TGapCosts, TTraceSwitch>
{
    typedef typename SubstituteAlignConfig<TAlignConfig>::Type TFreeEndGaps_;
    typedef AlignmentProfile<Global<TFreeEndGaps_>,  LinearGaps, TTraceSwitch> Type;
};

// profile for Gotoh
template <typename TAlignConfig, typename TGapCosts, typename TTraceSwitch>
struct SetUpAlignmentProfile<Gotoh, TAlignConfig, TGapCosts, TTraceSwitch>
{
    typedef typename SubstituteAlignConfig<TAlignConfig>::Type TFreeEndGaps_;
    typedef AlignmentProfile<Global<TFreeEndGaps_>, AffineGaps, TTraceSwitch> Type;
};

// profile for SmithWaterman
template <typename TAlignConfig, typename TGapCosts, typename TTraceSwitch>
struct SetUpAlignmentProfile<SmithWaterman, TAlignConfig, TGapCosts, TTraceSwitch>
{
    typedef typename SubstituteAlignConfig<TAlignConfig>::Type TFreeEndGaps_;
    typedef AlignmentProfile<Local<>,  TGapCosts, TTraceSwitch> Type;
};

// profile for WatermanEggert
template <typename TAlignConfig, typename TGapCosts, typename TTraceSwitch>
struct SetUpAlignmentProfile<WatermanEggert, TAlignConfig, TGapCosts, TTraceSwitch>
{
    typedef typename SubstituteAlignConfig<TAlignConfig>::Type TFreeEndGaps_;
    typedef AlignmentProfile<Local<WatermanEggert>, AffineGaps, TracebackSwitchedOn > Type;
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _doRunAlignment()
// ----------------------------------------------------------------------------

template <typename TTraceSegments, typename TSequenceHorizontal, typename TSequenceVertical, typename TScoreScheme,
        typename TBand, typename TAlignmentProfile>
inline typename Value<TScoreScheme>::Type
_doRunAlignment(TTraceSegments & traceSegments,
                TSequenceHorizontal const & seqHorizontal,
                TSequenceVertical const & seqVertical,
                TScoreScheme const & scoreScheme,
                TBand const & band,
                TAlignmentProfile const &  alignmentProfile)
{
    typedef typename GetGapSpec<TAlignmentProfile>::Type TGapSpec;
    typedef typename Value<TScoreScheme>::Type TScoreValue;
    typedef DPValue<TScoreValue, TGapSpec> TDPValue;
    typedef DPMatrix<TDPValue> TDPMatrix;

    typedef typename TraceBitMask::Type TTraceValue;
    typedef String<TTraceValue> TTraceMatrix;
    typedef typename Iterator<TTraceMatrix>::Type TTraceIterator;
    typedef typename IsTracebackOn<TAlignmentProfile>::Type TTracebackOn;
    typedef typename GetTrackerSpec<TAlignmentProfile>::Type TTrackerSpec;

    // construct trace Matrix
    TTraceMatrix traceMatrix;
    resize(traceMatrix, getTraceMatrixSize(seqHorizontal, seqVertical, band, TTracebackOn()), +TraceBitMask::NONE);

    // Traceback tracker
    Tracker<TScoreValue, TTraceIterator, TTrackerSpec> tracker;

    computeDPMatrix(tracker, traceMatrix, seqHorizontal, seqVertical, scoreScheme, band, alignmentProfile);

    computeTraceback(traceSegments, tracker, traceMatrix, seqHorizontal, seqVertical, band, alignmentProfile);

    return tracker._maxScore;
}


// ----------------------------------------------------------------------------
// Function _runAlignment()
// ----------------------------------------------------------------------------

template <typename TTraceSegments, typename TSequenceH, typename TSequenceV, typename TScoreScheme,
        typename TAlignConfig, typename TAlgoTag, typename TBand, typename TTraceSwitch>
inline typename Value<TScoreScheme>::Type
_runAlignment(TTraceSegments & traceSegments,
              TSequenceH const & seqHorizontal,
              TSequenceV const & seqVertical,
              TScoreScheme const & scoreScheme,
              TAlignConfig const & /*alignConfig*/,
              TAlgoTag const &  /*algoTag*/,
              TBand const & band,
              TTraceSwitch const & /*traceSwitch*/)
{
    SEQAN_ASSERT_GT(length(seqHorizontal), 0u);
    SEQAN_ASSERT_GT(length(seqVertical), 0u);

    if (scoreGapOpen(scoreScheme) == scoreGapExtend(scoreScheme))
    {
        typedef typename SetUpAlignmentProfile<TAlgoTag, TAlignConfig, LinearGaps, TTraceSwitch>::Type TAlignmentProfile;
        return _doRunAlignment(traceSegments, seqHorizontal, seqVertical, scoreScheme, band, TAlignmentProfile());
    }
    else
    {
        typedef typename SetUpAlignmentProfile<TAlgoTag, TAlignConfig, AffineGaps, TTraceSwitch>::Type TAlignmentProfile;
        return _doRunAlignment(traceSegments, seqHorizontal, seqVertical, scoreScheme, band, TAlignmentProfile());
    }
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_H_
