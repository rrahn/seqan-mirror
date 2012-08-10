/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
  ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ============================================================================
  Author: Rene Maerker <rene.maerker@fu-berlin.de>
  03.02.2011
  ============================================================================
  Implements the basic interface for storing the traceback of an alignment
  ==========================================================================*/

#ifndef STRING_SET_JOURNALED_JOIN_TRACE_H_
#define STRING_SET_JOURNALED_JOIN_TRACE_H_

namespace SEQAN_NAMESPACE_MAIN
{
	/**
	 * Encodes the directions of a trace back.
	 */
	enum TraceDirection
	{
			DIAGONAL,
			HORIZONTAL,
			VERTICAL
	};

	/**
	 * Basic interface for the JournalTraceDescriptor.
	 */
	template <typename TString>
	struct JournalTraceDescriptor;

	/**
	 * Specialization of the JournalTraceDescriptor for the use of Journal sequences.
	 * Describes the trace back in form of journal nodes in a sorted array.
	 * Note, that the nodes are entered in reversed order, because the trace back
	 * is parsed from the end to the beginning of the alignment.
	 */
	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	struct JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > >
	{
		typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString;
		typedef typename Position<TString>::Type TPos;
		typedef typename Size<TString>::Type TSize;
		typedef JournalEntry<TPos, TSize> TJournalNode;

		//stores the journal operations of the trace back reverse in consecutive order
		String<TJournalNode> revSortedOperation_;
		String<TValue, TBuffSpec> insertionBuffer_;

		JournalTraceDescriptor()
		{
			SEQAN_CHECKPOINT;
		}

		template <typename TPos>
		inline typename Reference<JournalTraceDescriptor>::Type
		operator[](TPos const pos)
		{
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<JournalTraceDescriptor const>::Type
		operator[](TPos const pos) const
		{
			return value(*this, pos);
		}

	};

	//_________________________META___________________________

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	struct Value<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >
	{
		typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
		typedef typename Position<TString_>::Type TPos_;
		typedef typename Size<TString_>::Type TSize_;
		typedef JournalEntry<TPos_, TSize_> Type;
	};

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	struct Value<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const>
	{
		typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
		typedef typename Position<TString_>::Type TPos_;
		typedef typename Size<TString_>::Type TSize_;
		typedef JournalEntry<TPos_, TSize_> Type;
	};

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	struct Reference<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >
	{
		typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
		typedef typename Value<JournalTraceDescriptor<TString_> >::Type & Type;
	};

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	struct Reference<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const>
	{
		typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
		typedef typename Value<JournalTraceDescriptor<TString_> >::Type const & Type;
	};

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	struct Size<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >
	{
		typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
		typedef typename Size<TString_>::Type Type;
	};

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	struct Size<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const>
		: Size<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >{};

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	struct Position<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >
	{
		typedef String<TValue, Journaled< THostSpec, TJournalSpec, TBuffSpec> > TString_;
		typedef typename Position<TString_>::Type Type;
	};

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	struct Position<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const>
		: Position<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > > {};


	//_____________________________Interface________________________________________

	/**
	 * Returns the trace nodes in order as they come.
	 */
	template <typename TString>
	inline String <typename Value<JournalTraceDescriptor<TString> >::Type> &
	getTrace(JournalTraceDescriptor<TString > & me)
	{
		SEQAN_CHECKPOINT;
		return me.revSortedOperation_;
	}
	/**
	 * Returns the trace nodes in order as they come.
	 */
	template <typename TString>
	inline String<typename Value<JournalTraceDescriptor<TString> >::Type> const &
	getTrace(JournalTraceDescriptor<TString> const & me)
	{
		SEQAN_CHECKPOINT;
		return me.revSortedOperation_;
	}

	/**
	 * Returns the trace in reverse (left-to-right) order. Use this function to access the nodes
	 * in consecutive order from the beginning to the end of the sequence.
	 */
	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	inline String <typename Value<JournalTraceDescriptor<String< TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >::Type>
	getTraceReverse(JournalTraceDescriptor<String< TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > & me)
	{
		SEQAN_CHECKPOINT;

		typedef String< TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournal;
		typedef typename Value<JournalTraceDescriptor<TJournal> >::Type  TEntry;
		typedef String<TEntry> TJournalEntries;

		TJournalEntries cpy(me.revSortedOperation_);
		reverse(cpy);
		return cpy;
	}

	/**
	 * Returns the trace in reverse (left-to-right) order. Use this function to access the nodes
	 * in consecutive order from the beginning to the end of the sequence.
	 */
	template <typename TString>
	inline String <typename Value<JournalTraceDescriptor<TString> const>::Type> const
	getTraceReverse(JournalTraceDescriptor<TString> const & me)
	{
		SEQAN_CHECKPOINT;
		return getTraceReverse(const_cast<JournalTraceDescriptor<TString> &>(me));
	}

	/**
	 * Returns the insertion buffer stored for the trace in case of insertions have been applied.
	 */
	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	inline String<TValue, TBuffSpec> &
	getInsertionBuffer(JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > & me)
	{
		SEQAN_CHECKPOINT;
		return me.insertionBuffer_;
	}

	/**
	 * Returns the insertion buffer stored for the trace in case of insertions have been applied.
	 */
	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	inline String<TValue, TBuffSpec> const &
	getInsertionBuffer(JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const & me)

	{
		SEQAN_CHECKPOINT;
		return me.insertionBuffer_;
	}

	/**
	 * Returns a reference to the node at the given position within the trace back.
	 */
	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	inline typename Reference<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >::Type
	value(JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > & me,
          typename Position<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >::Type const pos)
	{
		SEQAN_CHECKPOINT;
		return me.revSortedOperation_[pos];
	}

	/**
	 * Returns a reference to the node at the given position within the trace back.
	 */
	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	inline typename Reference<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const >::Type
	value(JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const & me,
		  typename Position<JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > >::Type const pos)
	{
		SEQAN_CHECKPOINT;
		return me.revSortedOperation_[pos];
	}

	/**
	 * Translates the current trace back operation into a journal operation that is stored as a journal entry
	 * in right-to-left order. Also the information of the horizontal operation are stored in an insertion buffer.
	 */
	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TStringSet, typename TId, typename TPos, typename TTraceValue>
	inline void
	_alignTracePrint(JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > & me,
			 TStringSet const & str,
			 TId const id1,
			 TPos const pos1,
			 TId const,
			 TPos const pos2,
			 TPos const segLen,
			 TTraceValue const tv)
	 {
		SEQAN_CHECKPOINT;

		typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TString;
		typedef typename Value<JournalTraceDescriptor<TString> >::Type TEntryString;
		typedef typename Value<TEntryString>::Type TJournalEntry;

		enum SegmentSource segmentSrc;
		TPos physicalPos;

		if (segLen == 0)
		{
			return;
		}
		switch (tv)
		{
		case DIAGONAL://matching area
			segmentSrc = SOURCE_ORIGINAL;
			physicalPos = pos2;
			break;
		case HORIZONTAL://insertion
			segmentSrc = SOURCE_PATCH;
			physicalPos = length(me.insertionBuffer_);

			append(me.insertionBuffer_, infix(str[id1],pos1, pos1 + segLen));
			break;
		case VERTICAL://deletion - nothing to be done here
			return;
			break;
		}
		appendValue(getTrace(me), TJournalEntry(segmentSrc, physicalPos, pos1, segLen));
	 }

	template <typename TStream, typename TJournal>
	TStream &
	operator <<(TStream & stream,
				JournalTraceDescriptor<TJournal> const & obj)
	{
		for (unsigned int i = 0; i < length(obj);++i)
		{
			if (obj[i].segmentSource == SOURCE_PATCH)
			{
				stream << obj[i] << " " <<infix(obj.insertionBuffer_,obj[i].physicalPosition, obj[i].physicalPosition + obj[i].length)<< std::endl;
			} else
			{
				stream << obj[i] << std::endl;
			}
		}
		return stream;
	}

	template <typename TJournal>
	inline typename Size<JournalTraceDescriptor<TJournal> >::Type
	length(JournalTraceDescriptor<TJournal> const & me)
	{
		SEQAN_CHECKPOINT;
		return length(me.revSortedOperation_);
	}

	/**
	 * Construct the tree balanced out of the balanced tree
	 */
	template <typename TValue, typename THostSpec, typename TBuffSpec>
	inline void
	constructAndSetJournalTree(String<TValue, Journaled<THostSpec, SortedArray, TBuffSpec> > & journalSeq,
							   JournalTraceDescriptor<String<TValue, Journaled<THostSpec, SortedArray, TBuffSpec> > > const & traceDescr)
	{
		SEQAN_CHECKPOINT;

		journalSeq._journalEntries._journalNodes = getTraceReverse(traceDescr);
		journalSeq._insertionBuffer = getInsertionBuffer(traceDescr);
	}

	template <typename TValue, typename THostSpec, typename TBuffSpec>
	inline void
	constructAndSetJournalTree(String<TValue, Journaled<THostSpec, UnbalancedTree, TBuffSpec> > & seq,
							   JournalTraceDescriptor<String<TValue, Journaled<THostSpec, UnbalancedTree, TBuffSpec> > > const & traceDescr)
	{
		SEQAN_CHECKPOINT;
		//implement tree construciton
//		seq._journalEntries = getTraceReverse(me);
//		seq._insertionBuffer = me.insertionBuffer_;
	}

}
#endif /* STRING_SET_JOURNALED_JOIN_TRACE_H_ */
