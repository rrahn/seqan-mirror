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

#ifndef CORE_INCLUDE_SEQAN_SEEDS2_SEEDS2_BANDED_CHAIN_ALIGNMENT_H_
#define CORE_INCLUDE_SEQAN_SEEDS2_BANDED_CHAIN_ALIGNMENT_H_

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

template <typename TGapSpec>
struct GetTrackerSpec<AlignmentProfile<Global<BandedChain>,  TGapSpec, TracebackSwitchedOn> >
{
    typedef BandedChain Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function bandedChainAlignment                                        [Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec, typename TSeedSpec, typename TSetSpec, typename TSetConfig,
        typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap, typename TBandSize>
typename TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps,
                     TBandSize const & minBandWidth)
{
    typedef typename Position<TSequence>::Type TPosition;
    typedef typename Size<TSequence>::Type TSize;
    typedef StringSet<String<TraceSegment<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;

    TScoreValue score = _runBandedChainAlignment(traceSet, source(row(align, 0)), source(row(align, 1)), seedSet,
                                             scoreSchemeAnchors, scoreSchemeGaps, minBandWidth);
    adapt(row(align,0), row(align,1), value(traceSet, 0));
    return score;
}

template <typename TSequence, typename TAlignSpec, typename TSeedSpec, typename TSetSpec, typename TSetConfig,
typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap>
typename TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps)
{
    return bandedChainAlignment(align, seedSet, scoreSchemeAnchors, scoreSchemeGaps, 15);
}

template <typename TSequence, typename TAlignSpec, typename TSeedSpec, typename TSetSpec, typename TSetConfig,
        typename TScoreValue, typename TScoreSpec, typename TBandSize>
typename TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     TBandSize const & minBandWidth)
{
    return bandedChainAlignment(align, seedSet, scoreScheme, scoreScheme, minBandWidth);
}

template <typename TSequence, typename TAlignSpec, typename TSeedSpec, typename TSetSpec, typename TSetConfig,
        typename TScoreValue, typename TScoreSpec>
typename TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return bandedChainAlignment(align, seedSet, scoreScheme, scoreScheme, 15);
}

// ----------------------------------------------------------------------------
// Function bandedChainAlignment                                         [Gaps]
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeedSpec,
        typename TSetSpec, typename TSetConfig,         typename TScoreValue, typename TScoreSpecAnchor,
        typename TScoreSpecGap, typename TBandSize>
typename TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps,
                     TBandSize const & minBandWidth)
{
    typedef typename Position<TSequenceH>::Type TPosition;
    typedef typename Size<TSequenceH>::Type TSize;
    typedef StringSet<String<TraceSegment<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;

    TScoreValue score = _runBandedChainAlignment(traceSet, source(gapsHorizontal), source(gapsVertical), seedSet,
                                             scoreSchemeAnchors, scoreSchemeGaps, minBandWidth);
    adapt(gapsHorizontal, gapsVertical, value(traceSet, 0));
    return score;
}

template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeedSpec,
        typename TSetSpec, typename TSetConfig,         typename TScoreValue, typename TScoreSpecAnchor,
        typename TScoreSpecGap>
typename TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps)
{
    return bandedChainAlignment(gapsHorizontal, gapsVertical, seedSet, scoreSchemeAnchors, scoreSchemeGaps, 15);
}

template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeedSpec,
        typename TSetSpec, typename TSetConfig, typename TScoreValue, typename TScoreSpec, typename TBandSize>
typename TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     TBandSize const & minBandWidth)
{
    return bandedChainAlignment(gapsHorizontal, gapsVertical, seedSet, scoreScheme, scoreScheme, minBandWidth);
}

template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeedSpec,
        typename TSetSpec, typename TSetConfig, typename TScoreValue, typename TScoreSpec>
typename TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return bandedChainAlignment(gapsHorizontal, gapsVertical, seedSet, scoreScheme, scoreScheme, 15);
}

// ----------------------------------------------------------------------------
// Function bandedChainAlignment                          [Graph<Alignment<> >]
// ----------------------------------------------------------------------------

template <typename TStringSet, typename TCargo, typename  TGraphSpec, typename TSeedSpec, typename TSetSpec,
       typename TSetConfig, typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap, typename TBandSize>
typename TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps,
                     TBandSize const & minBandWidth)
{
    typedef typename Position<TStringSet>::Type TPosition;
    typedef typename Size<TStringSet>::Type TSize;
    typedef StringSet<String<TraceSegment<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;

    TScoreValue score = _runBandedChainAlignment(traceSet, value(stringSet(alignmentGraph), 0),
                                                 value(stringSet(alignmentGraph), 1), seedSet, scoreSchemeAnchors,
                                                 scoreSchemeGaps, minBandWidth);

    adapt(alignmentGraph, positionToId(stringSet(alignmentGraph), 0), positionToId(stringSet(alignmentGraph), 1),
          value(traceSet, 0));

    return score;
}

template <typename TStringSet, typename TCargo, typename  TGraphSpec, typename TSeedSpec, typename TSetSpec,
       typename TSetConfig, typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap>
typename TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps)
{
    return bandedChainAlignment(alignmentGraph, seedSet, scoreSchemeAnchors, scoreSchemeGaps, 15);
}

template <typename TStringSet, typename TCargo, typename  TGraphSpec, typename TSeedSpec, typename TSetSpec,
       typename TSetConfig, typename TScoreValue, typename TScoreSpec, typename TBandSize>
typename TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     TBandSize const & minBandWidth)
{
    return bandedChainAlignment(alignmentGraph, seedSet, scoreScheme, scoreScheme, minBandWidth);
}

template <typename TStringSet, typename TCargo, typename  TGraphSpec, typename TSeedSpec, typename TSetSpec,
       typename TSetConfig, typename TScoreValue, typename TScoreSpec>
typename TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     SeedSet<TSeedSpec, TSetSpec, TSetConfig> const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme)
{
    return bandedChainAlignment(alignmentGraph, seedSet, scoreScheme, scoreScheme, 15);
}


}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_SEEDS2_BANDED_CHAIN_ALIGNMENT_H_
