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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_BASE_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_BASE_H_

namespace SEQAN_NAMESPACE_MAIN {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Global Alignment Algorithm Tags
// ----------------------------------------------------------------------------

// TODO(holtgrew): Rename MyersBitVector to Myers? Clashes with find module at the moment.
// TODO (rmaerker): which gap costs are used for Meyers and MeyersHirschberg
/**
.Tag.Pairwise Global Alignment Algorithms
..cat:Alignments
..summary:Tags used for selecting pairwise global alignment algorithms.
..tag
...Gotoh:Gotoh's for affine gap costs.
...NeedlemanWunsch:The Needleman-Wunsch algorithm for linear gap costs.
...Hirschberg:Hirschberg's algorithm using linear space.
...MyersBitVector:Myer's bit-vector algorithm.
...MyersHirschberg:Combination of Myer's and Hirschberg's algorithm.
..see:Function.globalAlignment
..see:Function.globalAlignmentScore
..include:seqan/align.h
*/

struct Gotoh_;
typedef Tag<Gotoh_> Gotoh;

struct NeedlemanWunsch_;
typedef Tag<NeedlemanWunsch_> NeedlemanWunsch;

struct Hirschberg_;
typedef Tag<Hirschberg_> Hirschberg;

struct MyersBitVector_;
typedef Tag<MyersBitVector_> MyersBitVector;

struct MyersHirschberg_;
typedef Tag<MyersHirschberg_> MyersHirschberg;

// ----------------------------------------------------------------------------
// Local Alignment Algorithm Tags
// ----------------------------------------------------------------------------

// TODO (rmaerker): which gap costs are used for local alignments
/**
.Tag.Pairwise Local Alignment Algorithms
..cat:Alignments
..summary:Tags used for selecting pairwise local alignment algorithms.
..tag
...SmithWaterman:Smith-Waterman algorithm for local alignments.
...WatermanEggert:Smith-Waterman algorithm with declumping to identify suboptimal local alignments.
..see:Function.localAlignment
..see:Class.LocalAlignmentEnumerator
..include:seqan/align.h
*/

struct SmithWaterman_;
typedef Tag<SmithWaterman_> SmithWaterman;

struct WatermanEggert_;
typedef Tag<WatermanEggert_> WatermanEggert;



// ----------------------------------------------------------------------------
// Gap Functions
// ----------------------------------------------------------------------------

/**
.Internal.Gap Costs
..cat:Alignments
..summary:Tags used for determine the used gap function.
..tag
...LinearGaps:Computes the alignment using linear gap costs.
...AffineGaps:Computes the alignment using affine gap costs.
..see:Tag.Global Alignment Algorithms
..see:Tag.Local Alignment Algorithms
..include:seqan/align.h
..remarks:For internal use only. See "Tag.Pairwise Global Alignment Algorithms"
and "Tag.Pairwise Local Alignment Algorithms" for algorithms that use different
gap costs. Use the scoring scheme to select different gap costs, but note that
NeddelmanWunsch and Gotoh overwrite the scoring scheme specialization.
*/

struct LinearGapFunction_;
typedef Tag<LinearGapFunction_> LinearGaps;

struct AffineGapFunction_;
typedef Tag<AffineGapFunction_> AffineGaps;


// ----------------------------------------------------------------------------
// Band
// ----------------------------------------------------------------------------

// The band can be either Off or On. If the band is wide one has to use
// On<Wide> instead of the default variant. This has to be determined at run time.

struct BandSwitchedOff_;
typedef Tag<BandSwitchedOff_> BandSwitchedOff;

struct WideBand_;
typedef Tag<WideBand_> WideBand;

template <typename TSpec = Default>
struct BandSwitchedOn;

template <typename TSpec = BandSwitchedOff>
struct Band{};

// ----------------------------------------------------------------------------
// Traceback
// ----------------------------------------------------------------------------

struct TracebackSwitchedOff_;
typedef Tag<TracebackSwitchedOff_> TracebackSwitchedOff;

struct TracebackSwitchedOn_;
typedef Tag<TracebackSwitchedOn_> TracebackSwitchedOn;

// traceback can be either On or Off.
template <typename TSpec = TracebackSwitchedOn>
struct Traceback
{};

// ----------------------------------------------------------------------------
// Traceback Encoding Table
// ----------------------------------------------------------------------------

// globally used bit mask to encode traceback directions
// Note, we only need three directions since we keep track of
// the origin of the optimal score per cell wihtin the dp value.

struct TraceBitMask
{
    typedef uint8_t Type;

    static const Type NONE = 0;                  //000000
    static const Type DIAGONAL = 1;              //100000
    static const Type VERTICAL = 2;              //010000
    static const Type HORIZONTAL = 4;            //001000
    static const Type HORIZONTAL_OPEN = 8;       //000100
    static const Type VERTICAL_OPEN = 16;        //000010
    static const Type FORBIDDEN = 32;            //000001
};

// ----------------------------------------------------------------------------
// DP Value
// ----------------------------------------------------------------------------

// The dp value is used to distinguish between linear and affine gap costs.
// It uses the configuration structures to forbid certain cells.
// This is used for the BandedChain and the WatermanEggert algorithm.

struct DPValueConfigDefault {};

struct DPValueConfigForbidden
{
    bool _forbidden;

    DPValueConfigForbidden() : _forbidden(false) {}

    DPValueConfigForbidden(DPValueConfigForbidden const & other) : _forbidden(other._forbidden){}

    DPValueConfigForbidden(bool const & forbidden) : _forbidden(forbidden){}
};

template <typename TScoreValue, typename TGapSpec, typename TConfig = DPValueConfigDefault>
struct DPValue{};

// ----------------------------------------------------------------------------
// Overlap Specialization
// ----------------------------------------------------------------------------

template <typename TFirstRowFree = False, typename TFirstColumnFree = False, typename TLastRowFree = False,
        typename TLastColumnFree = False>
struct FreeEndGaps{};

struct FirstRowAccessor_ {};
typedef Tag<FirstRowAccessor_> FirstRow;

struct FirstColumnAccessor_ {};
typedef Tag<FirstColumnAccessor_> FirstColumn;

struct LastRowAccessor_ {};
typedef Tag<LastRowAccessor_> LastRow;

struct LastColumnAccessor_ {};
typedef Tag<LastColumnAccessor_> LastColumn;


// ----------------------------------------------------------------------------
// Global algorithm
// ----------------------------------------------------------------------------

template <typename TSpec = FreeEndGaps<> >
struct Global{};

// ----------------------------------------------------------------------------
// Local Alignment Tag
// ----------------------------------------------------------------------------

template <typename TSpec = Default>
struct Local{};

// ----------------------------------------------------------------------------
// Split-Breakpoint algorithm
// ----------------------------------------------------------------------------

struct SplitBreakpoints_{};
typedef Tag<SplitBreakpoints_> SplitBreakpoint; // used for split breakpoint computation

// ----------------------------------------------------------------------------
// Banded-Chain algorithm
// ----------------------------------------------------------------------------

// TODO (rmaerker): move to seeds module.
struct BandedChain_{};
typedef Tag<BandedChain_> BandedChain; // used to compute banded chain alignment

// ----------------------------------------------------------------------------
// Alignment Profile
// ----------------------------------------------------------------------------

// The general alignment profile sets up the properties of the used dp algorithm
template <typename TAlignmentSpec, typename TGapsSpec = LinearGaps, typename TTracebackSpec = TracebackSwitchedOn >
struct AlignmentProfile{};

// ----------------------------------------------------------------------------
// DP Directions
// ----------------------------------------------------------------------------

// dp direction used for diagonal computation
struct DPDirectionDiagonal_;
typedef Tag<DPDirectionDiagonal_> DPDirectionDiagonal;

// dp direction used for vertical initialization
struct DPDirectionVertical_;
typedef Tag<DPDirectionVertical_> DPDirectionVertical;

// dp direction used for horizontal initialization
struct DPDirectionHorizontal_;
typedef Tag<DPDirectionHorizontal_> DPDirectionHorizontal;

// dp direction used in banded alignments if the first cell of a column is not
// at position 0;
struct DPDirectionUpperBand_;
typedef Tag<DPDirectionUpperBand_> DPDirectionUpperBand;

// dp direction used in banded alignments if the last cell of a column is not
// at position length(seqVertical);
struct DPDirectionLowerBand_;
typedef Tag<DPDirectionLowerBand_> DPDirectionLowerBand;

// standard dp direction while referring to all three neighboring cells
struct DPDirectionAll_;
typedef Tag<DPDirectionAll_> DPDirectionAll;

// used in local alignments or for initialization purposes
struct DPDirectionZero_;
typedef Tag<DPDirectionZero_> DPDirectionZero;

// special dp direction used for banded chain alignment. We simply adopt the
// score and the trace back direction of the current read value without
// computing it.
struct DPDirectionAdopt_;
typedef Tag<DPDirectionAdopt_> DPDirectionAdopt;


// ----------------------------------------------------------------------------
// DP Phases
// ----------------------------------------------------------------------------

// The dp algorithm runs column wise over dp matrix. Each column can be classified
// into one of the following four column types

// The first phase is the vertical initialization column, which is always the first
// column in a dp algorithm.
struct DPNoBandInit_;
typedef Tag<DPNoBandInit_> DPNoBandInitPhase; // Vertical Initialization

// In this phase the either the expansion of the band has reached the full column size
// or it is an unbanded dp algorithm. In both cases we have to distinguish between the
// dp direction of the banded and the unbanded dp algorithm within the first and the last
// cell of the column.
struct DPNoBandPhase_;
typedef Tag<DPNoBandPhase_> DPNoBandPhase; // The full column

// The first phase is the vertical initialization column, which is always the first
// column in a dp algorithm.
struct DPBandInit_;
typedef Tag<DPBandInit_> DPBandInitPhase; // Vertical Initialization

// In banded dp algorithms the second phase starts with a horizontal initialization.
// In this phase the full column size is approached while stepping forward in the
// dp matrix
struct DPBandFirstPhase_;
typedef Tag<DPBandFirstPhase_> DPBandFirstPhase;

// TODO(rmaerker): write me!
struct DPBandMiddlePhase_;
typedef Tag<DPBandMiddlePhase_> DPBandMiddlePhase;

// The last phase can only apply to banded dp algorithms, when the band reached the
// bottom of the dp matrix and begun to shrink.
struct DPBandLastPhase_;
typedef Tag<DPBandLastPhase_> DPBandLastPhase; // If Band the last part of the alignment

// The first phase is the vertical initialization column, which is always the first
// column in a dp algorithm.
struct DPWideBandInit_;
typedef Tag<DPWideBandInit_> DPWideBandInitPhase; // Vertical Initialization

// In banded dp algorithms the second phase starts with a horizontal initialization.
// In this phase the full column size is approached while stepping forward in the
// dp matrix
struct DPWideBandFirstPhase_;
typedef Tag<DPWideBandFirstPhase_> DPWideBandFirstPhase;

// TODO(rmaerker): write me!
struct DPWideBandMiddlePhase_;
typedef Tag<DPWideBandMiddlePhase_> DPWideBandMiddlePhase;

// The last phase can only apply to banded dp algorithms, when the band reached the
// bottom of the dp matrix and begun to shrink.
struct DPWideBandLastPhase_;
typedef Tag<DPWideBandLastPhase_> DPWideBandLastPhase; // If Band the last part of the alignment

// ----------------------------------------------------------------------------
// Cell Specifier
// ----------------------------------------------------------------------------

// Specifiers used to access the information of the different cell types per
// column.

struct FirstCell_;
typedef Tag<FirstCell_> FirstCell;

struct InnerCell_;
typedef Tag<InnerCell_> InnerCell;

struct LastCell_;
typedef Tag<LastCell_> LastCell;


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// GetAlingmentCore
// ----------------------------------------------------------------------------

template <typename TAlignmentProfile>
struct GetAlignmentSpec
{
    typedef Nothing Type;
};

template <typename TAlignmentSpec, typename TGapSpec, typename TTraceSpec>
struct GetAlignmentSpec<AlignmentProfile<TAlignmentSpec, TGapSpec, TTraceSpec> >
{
    typedef TAlignmentSpec Type;
};

template <typename TAlignmentSpec, typename TGapSpec, typename TTraceSpec>
struct GetAlignmentSpec<AlignmentProfile<TAlignmentSpec, TGapSpec, TTraceSpec> const >
    : GetAlignmentSpec<AlignmentProfile<TAlignmentSpec, TGapSpec, TTraceSpec> > {};

// ----------------------------------------------------------------------------
// GetGapSpec
// ----------------------------------------------------------------------------

template<typename TAlignmentProfile>
struct GetGapSpec
{
    typedef Nothing Type;
};

template<typename TAlignmentSpec, typename TGapSpec, typename TTraceSpec>
struct GetGapSpec< AlignmentProfile<TAlignmentSpec, TGapSpec, TTraceSpec> >
{
    typedef TGapSpec Type;
};

template<typename TAlignmentSpec, typename TGapSpec, typename TTraceSpec>
struct GetGapSpec< AlignmentProfile<TAlignmentSpec, TGapSpec, TTraceSpec> const >
        : GetGapSpec< AlignmentProfile<TAlignmentSpec, TGapSpec, TTraceSpec> >{};


// ----------------------------------------------------------------------------
// IsTracbackOn
// ----------------------------------------------------------------------------

template<typename TAlignmentProfile>
struct IsTracebackOn : True {};

template<typename TAlignmentSpec, typename TGapSpec>
struct IsTracebackOn<AlignmentProfile<TAlignmentSpec, TGapSpec, TracebackSwitchedOff > > : False {};

// ----------------------------------------------------------------------------
// IsFreeEndGap
// ----------------------------------------------------------------------------

template<typename TAlignmentSpec, typename TDPSide>
struct IsFreeEndGap : False{};

template<typename TAlignmentSpec, typename TGapSpec, typename TTracebackSpec, typename TDPSide>
struct IsFreeEndGap<AlignmentProfile< TAlignmentSpec, TGapSpec, TTracebackSpec>, TDPSide>
    : IsFreeEndGap<TAlignmentSpec, TDPSide>{};

template<typename TLocalSpec, typename TDPSide>
struct IsFreeEndGap<Local<TLocalSpec>, TDPSide> : True
{
};

template<typename TFirstColumn, typename TLastRow, typename TLastColumn>
struct IsFreeEndGap<Global<FreeEndGaps<True, TFirstColumn, TLastRow, TLastColumn> >, FirstRow> : True
{
};

template<typename TFirstRow, typename TLastRow, typename TLastColumn>
struct IsFreeEndGap<Global<FreeEndGaps<TFirstRow, True, TLastRow, TLastColumn> >, FirstColumn> : True
{
};

template<typename TFirstRow, typename TFirstColumn, typename TLastColumn>
struct IsFreeEndGap<Global<FreeEndGaps<TFirstRow, TFirstColumn, True, TLastColumn> >, LastRow> : True
{
};

template<typename TFirstRow, typename TFirstColumn, typename TLastRow>
struct IsFreeEndGap<Global<FreeEndGaps<TFirstRow, TFirstColumn, TLastRow, True> >, LastColumn> : True
{
};

// ----------------------------------------------------------------------------
// IsLocal
// ----------------------------------------------------------------------------

template <typename TAlignmentSpec>
struct IsLocal : False{};

template <typename TAlignmentSpec, typename TGapSpec, typename TTracebackSpec>
struct IsLocal<AlignmentProfile<TAlignmentSpec, TGapSpec, TTracebackSpec> >
    : IsLocal<TAlignmentSpec>{};


template <typename TLocalSpec>
struct IsLocal<Local< TLocalSpec > > : True{};


// ----------------------------------------------------------------------------
// IsGlobal
// ----------------------------------------------------------------------------

template <typename TAlignmentSpec>
struct IsGlobal : False{};

template <typename TAlignmentSpec, typename TGapSpec, typename TTracebackSpec>
struct IsGlobal<AlignmentProfile<TAlignmentSpec, TGapSpec, TTracebackSpec> >
: IsGlobal<TAlignmentSpec>{};


template <typename TGlobalSpec>
struct IsGlobal<Global< TGlobalSpec > > : True{};


// ----------------------------------------------------------------------------
// Spec
// ----------------------------------------------------------------------------

template <typename TAlignmentSpec, typename TGapSpec, typename TTracebackSpec>
struct Spec<AlignmentProfile<TAlignmentSpec, TGapSpec, TTracebackSpec> > : Spec<TAlignmentSpec> {};

template <typename TLocalSpec>
struct Spec<Local<TLocalSpec> >
{
    typedef TLocalSpec Type;
};

template <typename TGlobalSpec>
struct Spec<Global<TGlobalSpec> >
{
    typedef TGlobalSpec Type;
};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_BASE_H_
