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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): Switch to Host interface.

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAPS_BASE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAPS_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TSpec>
struct GapsIterator;

struct ArrayGaps_;
typedef Tag<ArrayGaps_> ArrayGaps;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Gaps
// ----------------------------------------------------------------------------

/**
.Class.Gaps
..implements:Concept.Sequence
..summary:Efficient storage of gaps for a sequence.
..signature:Gaps<TSequence, TSpec>
..description
...image:gaps_illustration|Illustration of Gaps object and positions with clipping.
..param.TSequence:The type of the underlying sequence.
...type:Concept.Sequence
..param.TSpec:Specialization tag.
..include:seqan/align.h
 */

template <typename TSequence, typename TSpec = ArrayGaps>
class Gaps;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

///.Metafunction.Value.param.T.type:Class.Gaps

template <typename TSequence, typename TSpec>
struct Value<Gaps<TSequence, TSpec> >
{
    typedef typename Value<TSequence>::Type           TAlphabet;
    typedef typename GappedValueType<TAlphabet>::Type Type;
};

template <typename TSequence, typename TSpec>
struct Value<Gaps<TSequence, TSpec> const> : Value<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

///.Metafunction.Iterator.param.T.type:Class.Gaps

template <typename TSequence, typename TSpec, typename TIteratorSpec>
struct Iterator<Gaps<TSequence, TSpec>, TIteratorSpec>
{
	typedef Iter<Gaps<TSequence, TSpec>, GapsIterator<TSpec> > Type;
};

template <typename TSequence, typename TSpec, typename TIteratorSpec>
struct Iterator<Gaps<TSequence, TSpec> const, TIteratorSpec>
{
	typedef Iter<Gaps<TSequence, TSpec> const, GapsIterator<TSpec> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

///.Metafunction.GetValue.param.T.type:Class.Gaps

template <typename TSequence, typename TSpec>
struct GetValue<Gaps<TSequence, TSpec> > : Value<Gaps<TSequence, TSpec> >
{};

template <typename TSequence, typename TSpec>
struct GetValue<Gaps<TSequence, TSpec> const> : GetValue<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

///.Metafunction.Position.param.T.type:Class.Gaps

template <typename TSequence, typename TSpec>
struct Position<Gaps<TSequence, TSpec> >
{
    typedef typename Position<TSequence>::Type TSeqPos_;
    typedef typename MakeSigned<TSeqPos_>::Type Type;
};

template <typename TSequence, typename TSpec>
struct Position<Gaps<TSequence, TSpec> const> : Position<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

///.Metafunction.Reference.param.T.type:Class.Gaps

template <typename TSequence, typename TSpec>
struct Reference<Gaps<TSequence, TSpec> >
{
	typedef typename Iterator<Gaps<TSequence, TSpec>, Standard>::Type TIterator_;
	typedef Proxy<IteratorProxy<TIterator_> > Type;
};

template <typename TSequence, typename TSpec>
struct Reference<Gaps<TSequence, TSpec> const>
{
	typedef typename Iterator<Gaps<TSequence, TSpec> const, Standard>::Type TIterator_;
	typedef Proxy<IteratorProxy<TIterator_> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

///.Metafunction.Size.param.T.type:Class.Gaps

template <typename TSequence, typename TSpec>
struct Size<Gaps<TSequence, TSpec> >
{
    typedef typename Size<TSequence>::Type Type;
};

template <typename TSequence, typename TSpec>
struct Size<Gaps<TSequence, TSpec> const> : Size<Gaps<TSequence, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Source
// ----------------------------------------------------------------------------

// TODO(holtgrew): Switch to Hosted Type interface

/**
.Metafunction.Source
..cat:Alignments
..summary:Return underlying sequence of Gaps/Alignments.
..signature:Source<T>::Type
..param.T:The type to query for underlying sequence.
..include:seqan/align.h
*/

///.Metafunction.Source.param.T.type:Class.Gaps

template <typename TSequence, typename TSpec>
struct Source<Gaps<TSequence, TSpec> >
{
    typedef TSequence Type;
};

template <typename TSequence, typename TSpec>
struct Source<Gaps<TSequence, TSpec> const> : Source<Gaps<TSequence, TSpec> >
{};

// TODO(holtgrew): Also prefix/suffix/infix? Should work!

// ----------------------------------------------------------------------------
// Metafunction IsSequence
// ----------------------------------------------------------------------------

///.Metafunction.IsSequence.param.T.type:Class.Gaps

template <typename TSequence, typename TSpec>
struct IsSequence<Gaps<TSequence, TSpec> >
{
    typedef True Type;
    static const bool VALUE = true;
};

template <typename TSequence, typename TSpec>
struct IsSequence<Gaps<TSequence, TSpec> const> : IsSequence<Gaps<TSequence, TSpec> >
{};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function toViewPosition()
// ----------------------------------------------------------------------------

/**
.Function.toViewPosition:
..summary:Transforms source to view position.
..cat:Alignments
..signature:toViewPosition(gaps, pos)
..param.gap:A Gaps object, e.g. a row in the alignment.
...type:Class.Gaps
..param.pos:Position in the original sequence to get the view position for.
..returns:The position in the view/gaps position.
..remarks:If $gap$ is a clipped alignment row, gaps in the clipped part will influence the result. The position $pos$ is counted from the unclipped begin position and must be greater or equal the clipped begin position of $gap$.
..see:Function.toSourcePosition
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function toSourcePosition()
// ----------------------------------------------------------------------------

/**
.Function.toSourcePosition:
..summary:Transforms view to source position, if the view position is a gap, the original position of the next non-gap entry is returned.
..cat:Alignments
..signature:toSourcePosition(gaps, pos)
..param.gap:A Gaps object, e.g. a row in the alignment.
...type:Class.Gaps
..param.pos:Position in the view sequence (this includes gaps) to get the original position for.
..returns:The position in the source sequence.
..remarks:If $gap$ is a clipped alignment row, gaps in the clipped part will influence the result. The position $pos$ is counted from the unclipped begin position.
..see:Function.toViewPosition
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function insertGap()
// ----------------------------------------------------------------------------

// Forward to removeGaps() which has to be implemented in each subclass.

template <typename TSequence, typename TSpec, typename TPosition>
inline void
insertGap(Gaps<TSequence, TSpec> & gaps, TPosition clippedViewPos)
{
    insertGaps(gaps, clippedViewPos, 1u);
}

// ----------------------------------------------------------------------------
// Function removeGap()
// ----------------------------------------------------------------------------

// Forward to removeGaps() which has to be implemented in each subclass.

template <typename TSequence, typename TSpec, typename TPosition>
inline typename Size<Gaps<TSequence, TSpec> >::Type
removeGap(Gaps<TSequence, TSpec> & gaps, TPosition clippedViewPos)
{
    return removeGaps(gaps, clippedViewPos, 1u);
}

// ----------------------------------------------------------------------------
// Function write()
// ----------------------------------------------------------------------------

template <typename TFile, typename TSource, typename TIDString, typename TSpec>
inline void
write(TFile & target,
	  Gaps<TSource, TSpec> const & source, 
	  TIDString const &,
	  Raw)
{
//IOREV _nodoc_ specialization not documented

	// Print gaps row
	typedef typename Iterator<Gaps<TSource, TSpec> const>::Type TIter;
	TIter begin_ = begin(source);
	TIter end_ = end(source);
	for (; begin_ != end_; ++begin_) {
		if (isGap(begin_))
			_streamPut(target, gapValue<char>());
		else 
			_streamPut(target, *begin_);
	}
}

// ----------------------------------------------------------------------------
// Function operator<<()                                      [stream operator]
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document appropriately.

template <typename TStream, typename TSource, typename TSpec>
inline TStream &
operator<<(TStream & stream, Gaps<TSource, TSpec> const & gaps)
{
    typedef Gaps<TSource, TSpec> const             TGaps;
    typedef typename Iterator<TGaps, Rooted>::Type TIter;
    typedef typename Value<TGaps>::Type            TGappedAlphabet;

    for (TIter it = begin(gaps, Rooted()); !atEnd(it); goNext(it))
    {
        // TODO(holtgrew): Ideally, we could simply print the expanded alphabet char but that is broken.
        if (isGap(it))
            stream << gapValue<char>();
        else
            stream << convert<char>(*it);
    }
    
    return stream;
}

// ----------------------------------------------------------------------------
// Function _pumpTraceToGaps()
// ----------------------------------------------------------------------------

// Internal function for converting AlignTrace<> objects into alignments in two Gaps objects.  Note that the traceback
// in the trace is stored in reverse, from back to front.  We insert the gaps in descending order of their position.
// The reason is that Gaps<> objects store the gaps in ascending order of coordinates in String<> objects and inserting
// at the end is in O(1) while inserting in the front is O(n).

template <typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV, typename TSize>
void _pumpTraceToGaps(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                      Gaps<TSequenceV, TGapsSpecV> & gapsV,
                      AlignTraceback<TSize> const & trace)
{
    typedef Gaps<TSequenceH, TGapsSpecH> TGapsH;
    typedef typename Iterator<TGapsH, Standard>::Type TGapsHIter;

    typedef Gaps<TSequenceV, TGapsSpecV> TGapsV;
    typedef typename Iterator<TGapsV, Standard>::Type TGapsVIter;

    // TODO(holtgrew): I don't understand the following.  Originally, this function used Align objects, but I did not understand it there either.
	// TODO(rausch): Pump trace into align_ (note: this is relatively slow code here. it could be improved if specialized to the Align Specs).
    clearGaps(gapsH);
    clearClipping(gapsH);
    clearGaps(gapsV);
    clearClipping(gapsV);

	TSize i = length(trace.sizes);  // Scan trace backwards.
	TGapsHIter itH = begin(gapsH);
	TGapsVIter itV = begin(gapsV);
	while (i > 0)
	{
		--i;
		TSize size = trace.sizes[i];
		switch ((int) trace.tvs[i])
		{
		case 1:  // Go horizontal.
			insertGaps(itV, size);
			break;

		case 2:  // Go vertical.
			insertGaps(itH, size);
			break;
		}
		goFurther(itH, size);
		goFurther(itV, size);
	}
}

// ----------------------------------------------------------------------------
// Function source()
// ----------------------------------------------------------------------------

/**
.Function.source
..summary:Return underlying sequence.
..cat:Alignments
..signature:source(obj)
..param.obj:The object to get underlying sequence of.
...type:Class.Gaps
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Function sourceSegment()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Rename/remove?

// We need some forwards for function sourceSegment().

template <typename TSequence, typename TSpec>
inline typename Position<Gaps<TSequence, TSpec> >::Type clippedBeginPosition(Gaps<TSequence, TSpec> const & gaps);
template <typename TSequence, typename TSpec>
inline typename Position<Gaps<TSequence, TSpec> >::Type clippedEndPosition(Gaps<TSequence, TSpec> const & gaps);
template <typename TSequence, typename TSpec, typename TPosition>
inline typename Position<TSequence>::Type
toSourcePosition(Gaps<TSequence, TSpec> const & gaps, TPosition clippedViewPos);

template <typename TSequence, typename TSpec>
inline typename Infix<TSequence>::Type
sourceSegment(Gaps<TSequence, TSpec> const & gaps)
{
    return infix(source(gaps), toSourcePosition(gaps, clippedBeginPosition(gaps)), toSourcePosition(gaps, clippedEndPosition(gaps)));
}

template <typename TSequence, typename TSpec>
inline typename Infix<TSequence>::Type
sourceSegment(Gaps<TSequence, TSpec> & gaps)
{
    return infix(source(gaps), toSourcePosition(gaps, clippedBeginPosition(gaps)), toSourcePosition(gaps, clippedEndPosition(gaps)));
}

// ----------------------------------------------------------------------------
// Function assignSource()
// ----------------------------------------------------------------------------

// TOOD(holtgrew): Switch to Hosted Type?

template <typename TSequence, typename TSpec, typename TValue>
inline void
assignSource(Gaps<TSequence, TSpec> & gaps, TValue const & value)
{
    assign(source(gaps), value);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
inline void clearGaps(Gaps<TSequence, TSpec> & gaps);
template <typename TSequence, typename TSpec>
inline void clearClipping(Gaps<TSequence, TSpec> & gaps);

template <typename TSequence, typename TSpec>
inline void clear(Gaps<TSequence, TSpec> & gaps)
{
    clearGaps(gaps);
    clearClipping(gaps);
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
inline bool operator==(Gaps<TSequence, TSpec> const & lhs,
                       Gaps<TSequence, TSpec> const & rhs)
{
    typename Comparator<Gaps<TSequence, TSpec> >::Type lex(lhs, rhs);
    return isEqual(lex);
}

template <typename TSequence, typename TSpec, typename TRightHandSide>
inline bool operator==(Gaps<TSequence, TSpec> const & lhs,
                       TRightHandSide const & rhs)
{
    typename Comparator<Gaps<TSequence, TSpec> >::Type lex(lhs, rhs);
    return isEqual(lex);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TSpec>
inline bool operator!=(Gaps<TSequence, TSpec> const & lhs,
                       Gaps<TSequence, TSpec> const & rhs)
{
    return !(lhs == rhs);
}


}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAPS_BASE_H_
