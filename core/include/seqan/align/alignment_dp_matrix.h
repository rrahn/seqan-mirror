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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_MATRIX_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_MATRIX_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(rmaerker): make dp matrix dependent from band.
template <typename TDPValue>
class DPMatrix
{
public:
    typedef String<TDPValue> TDataTable_;

    TDataTable_ _dataTable; //active recursion column which is used for computation

    DPMatrix() : _dataTable() {}

    template <typename TSequenceH, typename TSequenceV>
    DPMatrix(TSequenceH const & seqH, TSequenceV const & seqV) : _dataTable()
    {
        _initDPMatrix(*this, seqH, seqV, Band<BandSwitchedOff>());   //no band given so initialize without band
    }

    template <typename TSequenceH, typename TSequenceV, typename TBand>
    DPMatrix(TSequenceH const & seqH, TSequenceV const & seqV, TBand const & band) : _dataTable()
    {
        _initDPMatrix(*this, seqH, seqV, band);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Value
// ----------------------------------------------------------------------------

template <typename TDPValue>
struct Value<DPMatrix<TDPValue> >
{
    typedef TDPValue Type;
};

template <typename TDPValue>
struct Value<DPMatrix<TDPValue> const >
{
    typedef TDPValue const Type;
};

// ----------------------------------------------------------------------------
// Reference
// ----------------------------------------------------------------------------

template <typename TDPValue>
struct Reference<DPMatrix<TDPValue> >
{
    typedef TDPValue & Type;
};

template <typename TDPValue>
struct Reference<DPMatrix<TDPValue> const >
{
    typedef TDPValue const & Type;
};

// ----------------------------------------------------------------------------
// Size
// ----------------------------------------------------------------------------

template <typename TDPValue>
struct Size<DPMatrix<TDPValue> >
{
    typedef typename DPMatrix<TDPValue>::TDataTable_ TDataTable_;
    typedef typename Size<TDataTable_>::Type Type;
};

template <typename TDPValue>
struct Size<DPMatrix<TDPValue> const >
    : Size<DPMatrix<TDPValue> >{};

// ----------------------------------------------------------------------------
// Position
// ----------------------------------------------------------------------------

template <typename TDPValue>
struct Position<DPMatrix<TDPValue> >
{
    typedef typename DPMatrix<TDPValue>::TDataTable_ TDataTable_;
    typedef typename Position<TDataTable_>::Type Type;
};

template <typename TDPValue>
struct Position<DPMatrix<TDPValue> const >
    : Position<DPMatrix<TDPValue> >{};

// ----------------------------------------------------------------------------
// Iterator
// ----------------------------------------------------------------------------

template <typename TDPValue>
struct Iterator<DPMatrix<TDPValue>, Standard >
{
    typedef DPMatrix<TDPValue> TDPMatrix_;
    typedef typename TDPMatrix_::TDataTable_ TDataTable_;
    typedef typename Iterator<TDataTable_, Standard>::Type Type;
};

template <typename TDPValue>
struct Iterator<DPMatrix<TDPValue> const, Standard >
{
    typedef DPMatrix<TDPValue> TDPMatrix_;
    typedef typename TDPMatrix_::TDataTable_ TDataTable_;
    typedef typename Iterator<TDataTable_ const, Standard>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// function begin()
// ----------------------------------------------------------------------------

template <typename TDPValue>
inline typename Iterator<DPMatrix<TDPValue>, Standard>::Type
begin(DPMatrix<TDPValue> & recMatrix,
      Standard const & /*tag*/)
{
    return begin(recMatrix._dataTable);
}

template <typename TDPValue>
inline typename Iterator<DPMatrix<TDPValue> const, Standard>::Type
begin(DPMatrix<TDPValue> const & recMatrix,
      Standard const & /*tag*/)
{
    return begin(recMatrix._dataTable);
}

// ----------------------------------------------------------------------------
// function end()
// ----------------------------------------------------------------------------

template <typename TDPValue>
inline typename Iterator<DPMatrix<TDPValue>, Standard>::Type
end(DPMatrix<TDPValue> & recMatrix,
        Standard const & /*spec*/)
{
    return end(recMatrix._dataTable);
}

template <typename TDPValue>
inline typename Iterator<DPMatrix<TDPValue> const, Standard>::Type
end(DPMatrix<TDPValue> const & recMatrix,
        Standard const & /*spec*/)
{
    return end(recMatrix._dataTable);
}

// ----------------------------------------------------------------------------
// function _initDPMatrix()
// ----------------------------------------------------------------------------

template <typename TDPValue, typename TSequenceH, typename TSequenceV, typename TBand>
inline void _initDPMatrix(DPMatrix<TDPValue> & dpMatrix,
        TSequenceH const & seqH,
        TSequenceV const & seqV,
        TBand const & band)
{
    //single column matrix needs only one recursive column
    resize(dpMatrix._dataTable, getColumnSize(seqH, seqV, band), Exact());
}

//template <typename TDPValue, typename TSequenceH, typename TSequenceV, typename TBandSpec>
//inline void _initDPMatrix(DPMatrix<TDPValue> & dpMatrix,
//        TSequenceH const & seqH,
//        TSequenceV const & seqV,
//        Band<BandSwitchedOn<TBandSpec> > const & band)
//{
//    //single column matrix needs only one recursive column
//    resize(dpMatrix._dataTable, getColumnSize(seqH, seqV, band) + 1, Exact());
//}

// ----------------------------------------------------------------------------
// function _checkValidSequenceLengths()
// ----------------------------------------------------------------------------

template<typename TSeq0, typename TSeq1>
inline bool _checkValidSequenceLengths(TSeq0 const & seqH, TSeq1 const & seqV)
{
    return !empty(seqH) && !empty(seqV);
}

// ----------------------------------------------------------------------------
// function clear()
// ----------------------------------------------------------------------------

template <typename TDPValue>
inline void clear(DPMatrix<TDPValue> & dpMatrix)
{
    clear(dpMatrix._dataTable);
}

// ----------------------------------------------------------------------------
// function length()
// ----------------------------------------------------------------------------

template <typename TDPValue>
inline typename Size<DPMatrix<TDPValue> >::Type
length(DPMatrix<TDPValue> & dpMatrix)
{
    return length(dpMatrix._dataTable);
}

// ----------------------------------------------------------------------------
// function getColumnSize()
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TSequenceV, typename TSpec>
inline typename Size<TSequenceH>::Type
getColumnSize(TSequenceH const & /*seqH*/,
              TSequenceV const & seqV,
              Band<BandSwitchedOn<TSpec> > const & band)
{
    return _min(1 + static_cast<int>(length(seqV)), getBandSize(band));
}

template <typename TSequenceH, typename TSequenceV>
inline typename Size<TSequenceH>::Type
getColumnSize(TSequenceH const & /*seqH*/,
              TSequenceV const & seqV,
              Band<BandSwitchedOff> const & /*band*/)
{
    return length(seqV) + 1;
}

// ----------------------------------------------------------------------------
// function getRowSize()
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TSequenceV, typename TSpec>
inline typename Size<TSequenceH>::Type
getRowSize(TSequenceH const & seqH,
          TSequenceV const & /*seqV*/,
          Band<BandSwitchedOn<TSpec> > const & band)
{
    return _min(1 + static_cast<int>(length(seqH)), getBandSize(band));
}

template <typename TSequenceH, typename TSequenceV>
inline typename Size<TSequenceH>::Type
getRowSize(TSequenceH const & seqH,
              TSequenceV const & /*seqV*/,
              Band<BandSwitchedOff> const & /*band*/)
{
    return length(seqH) + 1;
}

// ----------------------------------------------------------------------------
// function getInitialColumnSize()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV>
inline typename Size<TSequenceH>::Type
getInitialColumnSize(TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Band<BandSwitchedOff> const & band)
{
    return getColumnSize(seqH, seqV, band);
}

template<typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline typename Size<TSequenceH>::Type
getInitialColumnSize(TSequenceH const & /*seqH*/,
                      TSequenceV const & seqV,
                      Band<BandSwitchedOn<TBandSpec> > const & band)
{
    // The minimal length of a column is one.
    return _max(1, _min(0, getUpperDiagonal(band)) + _min(1 + static_cast<int>(length(seqV)), 1 - getLowerDiagonal(band)));
}

// ----------------------------------------------------------------------------
// function getInitialRowSize()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV>
inline typename Size<TSequenceH>::Type
getInitialRowSize(TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Band<BandSwitchedOff> const & band)
{
    return getRowSize(seqH, seqV, band);
}

template<typename TSequenceH, typename TSequenceV, typename TBandSpec>
inline typename Size<TSequenceH>::Type
getInitialRowSize(TSequenceH const & seqH,
                      TSequenceV const & /*seqV*/,
                      Band<BandSwitchedOn<TBandSpec> > const & band)
{
    // The minimal length of a row is one.
    return _max(1,_min(1 + static_cast<int>(length(seqH)), 1 + getUpperDiagonal(band)));
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_DP_MATRIX_H_
