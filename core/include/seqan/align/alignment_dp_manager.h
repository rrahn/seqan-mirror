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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_MANAGER_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_MANAGER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TAlignmentProfile, typename TDPPhase>
struct SetUpColumnManager
{
    typedef Nothing Type;
};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TDPDirection, typename TTrackScore = False>
struct Cell {};

template <typename TFirstCell_, typename TInnerCell_, typename TLastCell_>
struct ColumnManager
{
    typedef TFirstCell_ TFirstCell;
    typedef TInnerCell_ TInnerCell;
    typedef TLastCell_ TLastCell;
};


template <typename TAlignmentProfile, typename TBand>
class DPManager
{
public:

    int _spanDp;
    int _spanTrace;
    int _spanSeqVBegin;
    int _spanSeqVEnd;

    DPManager() : _spanDp(0),
                      _spanTrace(0),
                      _spanSeqVBegin(0),
                      _spanSeqVEnd(0){}

    DPManager(DPManager const & other) : _spanDp(other._spanDp),
                                                 _spanTrace(other._spanTrace),
                                                 _spanSeqVBegin(other._spanSeqVBegin),
                                                 _spanSeqVEnd(other._spanSeqVEnd){}

    template <typename TSize>
    DPManager(TSize const & columnSize, TSize const & initialColumnSize) : _spanDp(),
                    _spanTrace(),
                    _spanSeqVBegin(),
                    _spanSeqVEnd()
    {
        _init(*this, columnSize, initialColumnSize);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

//template <typename TDPManager, typename TColType>
//struct GetColumnManager
//{
//    typedef Nothing Type;
//};
//
//template <typename TDPManager>
//struct GetColumnManager<TDPManager, InitialColumn>
//{
//    typedef typename TDPManager::TInitColumnManager_ Type;
//};
//
//template <typename TDPManager>
//struct GetColumnManager<TDPManager, BandOpenColumn>
//{
//    typedef typename TDPManager::TBeginColumnManager_ Type;
//};
//
//template <typename TDPManager>
//struct GetColumnManager<TDPManager, FullColumn>
//{
//    typedef typename TDPManager::TMiddleColumnManager_ Type;
//};
//
//template <typename TDPManager>
//struct GetColumnManager<TDPManager, BandCloseColumn>
//{
//    typedef typename TDPManager::TEndColumnManager_ Type;
//};

// TODO (rmaerker): move to banded chain algorithm implementation
// cell manager for the banded chain alignment

template <typename TGapCostSpec, typename TTracebackSpec>
struct SetUpColumnManager<AlignmentProfile<Global<BandedChain>, TGapCostSpec, TTracebackSpec>,
        DPNoBandInitPhase>
{
        // defines the specialization of the first cell
    typedef Cell<DPDirectionAdopt, False> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAdopt, False> TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionAdopt, True> TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TGapCostSpec, typename TTracebackSpec>
struct SetUpColumnManager<AlignmentProfile<Global<BandedChain>, TGapCostSpec, TTracebackSpec>, DPNoBandPhase>
{
        // defines the specialization of the first cell
    typedef Cell<DPDirectionAdopt, False> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, False> TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionAll, True> TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TGapCostSpec, typename TTracebackSpec>
struct SetUpColumnManager<AlignmentProfile<Global<BandedChain>, TGapCostSpec, TTracebackSpec>, DPBandInitPhase>
{
    // defines the specialization of the first cell
    typedef Cell<DPDirectionAdopt, False> TFirstCell_;

    // defines the specialization of the inner cells
    typedef Cell<DPDirectionAdopt, False> TInnerCell_;

    // defines the specialization of the last cell
    typedef Cell<DPDirectionAdopt, False> TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TGapCostSpec, typename TTracebackSpec>
struct SetUpColumnManager<AlignmentProfile<Global<BandedChain>, TGapCostSpec, TTracebackSpec>, DPBandFirstPhase>
{
    // defines the specialization of the first cell
    typedef Cell<DPDirectionAdopt, False> TFirstCell_;

    // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, False> TInnerCell_;

    // defines the specialization of the last cell
    typedef Cell<DPDirectionLowerBand, False> TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TGapCostSpec, typename TTracebackSpec>
struct SetUpColumnManager<AlignmentProfile<Global<BandedChain>, TGapCostSpec, TTracebackSpec>, DPBandMiddlePhase>
{
        // defines the specialization of the first cell
    typedef Cell<DPDirectionUpperBand, False> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, False> TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionLowerBand, False> TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TGapCostSpec, typename TTracebackSpec>
struct SetUpColumnManager<AlignmentProfile<Global<BandedChain>, TGapCostSpec, TTracebackSpec>, DPBandLastPhase>
{
        // defines the specialization of the first cell
    typedef Cell<DPDirectionUpperBand, False> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, False> TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionAll, True> TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TGapCostSpec, typename TTracebackSpec>
struct SetUpColumnManager<AlignmentProfile<Global<BandedChain>, TGapCostSpec, TTracebackSpec>, DPWideBandInitPhase>
{
    // defines the specialization of the first cell
    typedef Cell<DPDirectionAdopt, False> TFirstCell_;

    // defines the specialization of the inner cells
    typedef Cell<DPDirectionAdopt, False> TInnerCell_;

    // defines the specialization of the last cell
    typedef Cell<DPDirectionAdopt, False> TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TGapCostSpec, typename TTracebackSpec>
struct SetUpColumnManager<AlignmentProfile<Global<BandedChain>, TGapCostSpec, TTracebackSpec>, DPWideBandFirstPhase>
{
    // defines the specialization of the first cell
    typedef Cell<DPDirectionAdopt, False> TFirstCell_;

    // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, False> TInnerCell_;

    // defines the specialization of the last cell
    typedef Cell<DPDirectionLowerBand, False> TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TGapCostSpec, typename TTracebackSpec>
struct SetUpColumnManager<AlignmentProfile<Global<BandedChain>, TGapCostSpec, TTracebackSpec>, DPWideBandMiddlePhase>
{
        // defines the specialization of the first cell
    typedef Cell<DPDirectionAdopt, False> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, False> TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionAll, True> TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TGapCostSpec, typename TTracebackSpec>
struct SetUpColumnManager<AlignmentProfile<Global<BandedChain>, TGapCostSpec, TTracebackSpec>, DPWideBandLastPhase>
{
        // defines the specialization of the first cell
    typedef Cell<DPDirectionUpperBand, False> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, False> TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionAll, True> TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};


// setup column manager for standard dp algorithms

// construction of cell manager for vertical initialization
template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPNoBandInitPhase>
{
    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;
    typedef typename IsFreeEndGap<TAlignmentProfile, FirstColumn>::Type TFirstColumn_;
    typedef typename Or<IsLocal<TAlignmentProfile>,
        IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

    // defines the specialization of the first cell
    typedef Cell<DPDirectionZero, TIsTrackable_> TFirstCell_;

    // defines the specialization of the inner cells
    typedef Cell<typename If<Or<IsLocal<TAlignmentProfile>, TFirstColumn_>, DPDirectionZero, DPDirectionVertical>::Type,
                TIsTrackable_ > TInnerCell_;

    // defines the specialization of the last cell
    typedef Cell<typename If<Or<IsLocal<TAlignmentProfile>, TFirstColumn_>, DPDirectionZero, DPDirectionVertical>::Type,
    typename Or<TIsTrackable_, IsFreeEndGap<TAlignmentProfile, LastRow> >::Type >  TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPNoBandPhase>
{
    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;
    typedef typename IsFreeEndGap<TAlignmentProfile, FirstRow>::Type TFirstRow_;

    typedef typename Or<IsLocal<TAlignmentProfile>,
                IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint>  >::Type TIsTrackable_;

        // defines the specialization of the first cell
    typedef Cell<typename If<Or<IsLocal<TAlignmentProfile>, TFirstRow_>, DPDirectionZero, DPDirectionHorizontal>::Type, TIsTrackable_> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, TIsTrackable_ > TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionAll, typename Or<TIsTrackable_, IsFreeEndGap<TAlignmentProfile, LastRow> >::Type > TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell for initialization
};

// The column manager used for normal sized bands.

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPBandInitPhase>
{
    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;
    typedef typename IsFreeEndGap<TAlignmentProfile, FirstColumn>::Type TFirstColumn_;
    typedef typename Or<IsLocal<TAlignmentProfile>,
    IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

        // defines the specialization of the first cell
    typedef Cell<DPDirectionZero, TIsTrackable_> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<typename If<Or<IsLocal<TAlignmentProfile>, TFirstColumn_>, DPDirectionZero, DPDirectionVertical>::Type,
    TIsTrackable_ > TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<typename If<Or<IsLocal<TAlignmentProfile>, TFirstColumn_>, DPDirectionZero, DPDirectionVertical>::Type,
                 TIsTrackable_>  TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};


template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPBandFirstPhase>
{

    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;
    typedef typename IsFreeEndGap<TAlignmentProfile, FirstRow>::Type TFirstRow_;

    typedef typename Or<IsLocal<TAlignmentProfile>,
            IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

    // defines the specialization of the first cell
    typedef Cell<typename If<Or<IsLocal<TAlignmentProfile>, TFirstRow_>, DPDirectionZero, DPDirectionHorizontal>::Type, TIsTrackable_> TFirstCell_;

    // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, TIsTrackable_ > TInnerCell_;

    // defines the specialization of the last cell
    typedef Cell<DPDirectionLowerBand, TIsTrackable_ > TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell for initialization
};

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPBandMiddlePhase>
{

    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;

    typedef typename Or<IsLocal<TAlignmentProfile>,
                    IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

        // defines the specialization of the first cell
    typedef Cell<DPDirectionUpperBand, TIsTrackable_> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, TIsTrackable_ > TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionLowerBand, TIsTrackable_ > TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell for initialization
};

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPBandLastPhase>
{
    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;

    typedef typename Or<IsLocal<TAlignmentProfile>,
        IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

        // defines the specialization of the first cell
    typedef Cell<DPDirectionUpperBand, TIsTrackable_> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, TIsTrackable_ > TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionAll, typename Or<TIsTrackable_, IsFreeEndGap<TAlignmentProfile, LastRow> >::Type > TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell for initialization
};

// The column manager used for wide bands.

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPWideBandInitPhase>
{
    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;
    typedef typename IsFreeEndGap<TAlignmentProfile, FirstColumn>::Type TFirstColumn_;
    typedef typename Or<IsLocal<TAlignmentProfile>,
    IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

        // defines the specialization of the first cell
    typedef Cell<DPDirectionZero, TIsTrackable_> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<typename If<Or<IsLocal<TAlignmentProfile>, TFirstColumn_>, DPDirectionZero, DPDirectionVertical>::Type,
    TIsTrackable_ > TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<typename If<Or<IsLocal<TAlignmentProfile>, TFirstColumn_>, DPDirectionZero, DPDirectionVertical>::Type,
                 TIsTrackable_>  TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell manager for initialization
};

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPWideBandFirstPhase>
{

    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;
    typedef typename IsFreeEndGap<TAlignmentProfile, FirstRow>::Type TFirstRow_;

    typedef typename Or<IsLocal<TAlignmentProfile>,
            IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

    // defines the specialization of the first cell
    typedef Cell<typename If<Or<IsLocal<TAlignmentProfile>, TFirstRow_>, DPDirectionZero, DPDirectionHorizontal>::Type, TIsTrackable_> TFirstCell_;

    // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, TIsTrackable_ > TInnerCell_;

    // defines the specialization of the last cell
    typedef Cell<DPDirectionLowerBand, TIsTrackable_ > TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell for initialization
};

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPWideBandMiddlePhase>
{
    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;
    typedef typename IsFreeEndGap<TAlignmentProfile, FirstRow>::Type TFirstRow_;

    typedef typename Or<IsLocal<TAlignmentProfile>,
                IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint>  >::Type TIsTrackable_;

        // defines the specialization of the first cell
    typedef Cell<typename If<Or<IsLocal<TAlignmentProfile>, TFirstRow_>, DPDirectionZero, DPDirectionHorizontal>::Type, TIsTrackable_> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, TIsTrackable_ > TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionAll, typename Or<TIsTrackable_, IsFreeEndGap<TAlignmentProfile, LastRow> >::Type > TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell for initialization
};

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPWideBandLastPhase>
{
    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;

    typedef typename Or<IsLocal<TAlignmentProfile>,
        IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

        // defines the specialization of the first cell
    typedef Cell<DPDirectionUpperBand, TIsTrackable_> TFirstCell_;

        // defines the specialization of the inner cells
    typedef Cell<DPDirectionAll, TIsTrackable_ > TInnerCell_;

        // defines the specialization of the last cell
    typedef Cell<DPDirectionAll, typename Or<TIsTrackable_, IsFreeEndGap<TAlignmentProfile, LastRow> >::Type > TLastCell_;

    typedef ColumnManager<TFirstCell_, TInnerCell_, TLastCell_> Type; // constructs the cell for initialization
};

// The column manager used for one-size bands.

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPSmallBandInit>
{
    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;
    typedef typename Or<IsLocal<TAlignmentProfile>,
    IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

     // defines the specialization of the first cell
    typedef Cell<DPDirectionZero, TIsTrackable_> TCell_;

    // Small bands have only a width of one, so there is only one cell per column.
    typedef ColumnManager<TCell_, TCell_, TCell_> Type;
};

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPSmallBandDetached>
{

    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;
    typedef typename IsFreeEndGap<TAlignmentProfile, FirstRow>::Type TFirstRow_;

    typedef typename Or<IsLocal<TAlignmentProfile>,
            IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

    typedef Cell<DPDirectionDiagonal, TIsTrackable_> TCell_;
    // Small bands have only a width of one, so there is only one cell per column.
    typedef ColumnManager<TCell_, TCell_, TCell_> Type;
};

template <typename TAlignmentProfile>
struct SetUpColumnManager<TAlignmentProfile, DPSmallBandEnd>
{
    typedef typename GetAlignmentSpec<TAlignmentProfile>::Type TAlignmentSpec_;

    typedef typename Or<IsLocal<TAlignmentProfile>,
        IsSameType<typename Spec<TAlignmentSpec_>::Type, SplitBreakpoint> >::Type TIsTrackable_;

        // defines the specialization of the first cell
    typedef Cell<DPDirectionDiagonal,
            typename Or<TIsTrackable_, IsFreeEndGap<TAlignmentProfile, LastRow> >::Type> TCell_;

    typedef ColumnManager<TCell_, TCell_, TCell_> Type; // constructs the cell for initialization
};

// ----------------------------------------------------------------------------
// Metafunction GetDPDirection
// ----------------------------------------------------------------------------

template <typename TCell>
struct GetDpDirection{};

template <typename TDPDirectionSpec, typename TTraceSpec>
struct GetDpDirection<Cell<TDPDirectionSpec, TTraceSpec> >
{
    typedef TDPDirectionSpec Type;
};

template <typename TDPDirectionSpec, typename TTraceSpec>
struct GetDpDirection<Cell<TDPDirectionSpec, TTraceSpec> const >
{
    typedef TDPDirectionSpec const Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsToTrack
// ----------------------------------------------------------------------------

template <typename TCell>
struct IsToTrack : False {};

template <typename TDPDirectionSpec>
struct IsToTrack<Cell<TDPDirectionSpec, True> > : True {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _init()                                             [band disabled]
// ----------------------------------------------------------------------------

// initializes the column director given the column size and the initial column size
template <typename TAlignmentProfile, typename TSize>
inline void
_init(DPManager<TAlignmentProfile, Band<BandSwitchedOff> > & columnManager,
      TSize const & columnSize,
      TSize const & /*initialColumnSize*/)
{
    columnManager._spanDp = 1 - columnSize; // move to begin of column
    columnManager._spanTrace = 1; // move always one step forward
    columnManager._spanSeqVBegin = 0;
    columnManager._spanSeqVEnd = columnSize -2;
}

// ----------------------------------------------------------------------------
// Function _init()                                              [band enabled]
// ----------------------------------------------------------------------------

template <typename TAlignmentProfile, typename TBandSpec, typename TSize>
inline void
_init(DPManager<TAlignmentProfile, Band<BandSwitchedOn<TBandSpec> > > & columnManager,
      TSize const & columnSize,
      TSize const & initialColumnSize)
{
    columnManager._spanDp = 1 - initialColumnSize;
    columnManager._spanTrace = 1 + columnSize - initialColumnSize;
    columnManager._spanSeqVBegin = -1;
    columnManager._spanSeqVEnd = initialColumnSize -2;
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TMatrix, typename TAlignProfile, typename TBand>
inline typename Iterator<TMatrix>::Type
begin(TMatrix & matrix, DPManager<TAlignProfile, TBand> const & colManager)
{
    return begin(matrix, colManager, Standard());
}

template <typename TMatrix, typename TAlignProfile, typename TBand>
inline typename Iterator<TMatrix>::Type
begin(TMatrix const & matrix, DPManager<TAlignProfile, TBand> const & colManager)
{
    return begin(matrix, colManager, Standard());
}

template <typename TDPValue, typename TAlignProfile, typename TBand, typename TIterSpec>
inline typename Iterator<DPMatrix<TDPValue> const >::Type
begin(DPMatrix<TDPValue> const & dpMatrix,
      DPManager<TAlignProfile, TBand> & colManager,
      TIterSpec const & /*spec*/)
{
    return end(dpMatrix, TIterSpec()) + (colManager._spanDp - 1);
}

template <typename TDPValue, typename TAlignProfile, typename TBand, typename TIterSpec>
inline typename Iterator<DPMatrix<TDPValue> const >::Type
begin(DPMatrix<TDPValue> const & dpMatrix,
      DPManager<TAlignProfile, TBand > const & colManager,
      TIterSpec const & /*spec*/)
{
    return end(dpMatrix, TIterSpec()) + (colManager._spanDp - 1);
}

//template <typename TDPValue, typename TAlignProfile, typename TBandSpec, typename TIterSpec>
//inline typename Iterator<DPMatrix<TDPValue> const >::Type
//begin(DPMatrix<TDPValue> const & dpMatrix,
//      DPManager<TAlignProfile, Band<BandSwitchedOn<TBandSpec> > > & colManager,
//      TIterSpec const & /*spec*/)
//{
//    return end(dpMatrix, TIterSpec()) + (colManager._spanDp - 1);
//}
//
//template <typename TDPValue, typename TAlignProfile, typename TBandSpec, typename TIterSpec>
//inline typename Iterator<DPMatrix<TDPValue> const >::Type
//begin(DPMatrix<TDPValue> const & dpMatrix,
//      DPManager<TAlignProfile, Band<BandSwitchedOn<TBandSpec> > > const & colManager,
//      TIterSpec const & /*spec*/)
//{
//    return end(dpMatrix, TIterSpec()) + (colManager._spanDp - 1);
//}

template <typename TTraceMatrix, typename TDPManager, typename TIterSpec>
inline typename Iterator<TTraceMatrix>::Type
begin(TTraceMatrix & traceMatrix, TDPManager const & colManager, TIterSpec const & /*spec*/)
{
    return begin(traceMatrix, TIterSpec()) + (colManager._spanTrace - 1);
}

template <typename TTraceMatrix, typename TDPManager, typename TIterSpec>
inline typename Iterator<TTraceMatrix const>::Type
begin(TTraceMatrix const & traceMatrix, TDPManager const & colManager, TIterSpec const & /*spec*/)
{
    return begin(traceMatrix, TIterSpec()) + (colManager._spanTrace - 1);
}

// ----------------------------------------------------------------------------
// Function _correctSpan()
// ----------------------------------------------------------------------------

// In the banded version we need to correct the span once when we go to the last
// DP phase to refer to the correct cells in the next column.

template <typename TAlignmentProfile, typename TBand>
inline void
_correctSpan(DPManager<TAlignmentProfile, TBand> const & /*dpManager*/)
{
    // nothing to do
}

template <typename TAlignmentProfile, typename TBandSpec>
inline void
_correctSpan(DPManager<TAlignmentProfile, Band<BandSwitchedOn<TBandSpec> > > & dpManager)
{
   --dpManager._spanDp;
   --dpManager._spanTrace;
}

//template <typename TAlignmentProfile, typename TBand>
//inline void
//_correctSpan(DPManager<TAlignmentProfile, TBand> & columnManager,
//             DPBandLastPhase const & /*colType*/)
//{
//   --columnManager._spanDp;
//   --columnManager._spanTrace;
//}
//
//template <typename TAlignmentProfile, typename TBandSpec>
//inline void
//_correctSpan(DPManager<TAlignmentProfile, Band<BandSwitchedOn<TBandSpec> > > & columnManager,
//             DPWideBandMiddlePhase const & /*colType*/)
//{
//   --columnManager._spanDp;
//   --columnManager._spanTrace;
//}

// ----------------------------------------------------------------------------
// Function update()                                                  [no band]
// ----------------------------------------------------------------------------

template <typename TAlignmentProfile, typename TBand, typename TColumnType>
inline void
update(DPManager<TAlignmentProfile, TBand> const & /*columnManager*/, TColumnType const & /*colType*/)
{
    // nothing to do here
}

// ----------------------------------------------------------------------------
// Function update                                                [normal band]
// ----------------------------------------------------------------------------

template <typename TAlignmentProfile, typename TBandSpec>
inline void
update(DPManager<TAlignmentProfile, Band<BandSwitchedOn<TBandSpec> > > & columnManager,
       DPBandFirstPhase const & /*colType*/)
{
    ++columnManager._spanSeqVEnd;
    --columnManager._spanDp;
    --columnManager._spanTrace;
}

template <typename TAlignmentProfile, typename TBandSpec>
inline void
update(DPManager<TAlignmentProfile, Band<BandSwitchedOn<TBandSpec> > > & columnManager,
       DPBandMiddlePhase const & /*colType*/)
{
    ++columnManager._spanSeqVBegin;
    ++columnManager._spanSeqVEnd;
}

template <typename TAlignmentProfile, typename TBandSpec>
inline void
update(DPManager<TAlignmentProfile, Band<BandSwitchedOn<TBandSpec> > > & columnManager,
       DPBandLastPhase const & /*colType*/)
{
    ++columnManager._spanSeqVBegin;
    ++columnManager._spanDp;
    ++columnManager._spanTrace;
}

// ----------------------------------------------------------------------------
// Function update                                                  [wide band]
// ----------------------------------------------------------------------------

template <typename TAlignmentProfile, typename TBandSpec>
inline void
update(DPManager<TAlignmentProfile, Band<BandSwitchedOn<TBandSpec> > > & columnManager,
       DPWideBandFirstPhase const & /*colType*/)
{
    ++columnManager._spanSeqVEnd;
    --columnManager._spanDp;
    --columnManager._spanTrace;
}

template <typename TAlignmentProfile, typename TBandSpec>
inline void
update(DPManager<TAlignmentProfile, Band<BandSwitchedOn<TBandSpec> > > & columnManager,
       DPWideBandLastPhase const & /*colType*/)
{
    ++columnManager._spanSeqVBegin;
    ++columnManager._spanDp;
    ++columnManager._spanTrace;
}

// ----------------------------------------------------------------------------
// Function update                                                 [small band]
// ----------------------------------------------------------------------------

template <typename TAlignmentProfile, typename TBandSpec>
inline void
update(DPManager<TAlignmentProfile, Band<BandSwitchedOn<TBandSpec> > > & columnManager,
       DPSmallBandDetached const & /*colType*/)
{
    ++columnManager._spanSeqVBegin;
    ++columnManager._spanSeqVEnd;
}

template <typename TAlignmentProfile, typename TBandSpec>
inline void
update(DPManager<TAlignmentProfile, Band<BandSwitchedOn<TBandSpec> > > & columnManager,
       DPSmallBandEnd const & /*colType*/)
{
    ++columnManager._spanSeqVBegin;
    ++columnManager._spanSeqVEnd;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_MANAGER_H_
