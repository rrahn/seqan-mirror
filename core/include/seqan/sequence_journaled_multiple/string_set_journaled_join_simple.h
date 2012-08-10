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
  02.02.2011
  ============================================================================
  todo enter description of this file
  ==========================================================================*/

#ifndef STRING_SET_JOURNALED_JOIN_SIMPLE_H_
#define STRING_SET_JOURNALED_JOIN_SIMPLE_H_

namespace SEQAN_NAMESPACE_MAIN
{



	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TScore>
	void
	internalJoin_(String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journal,
				  String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & refJournal,
				  TScore const &,
				  SimpleJoin const &)
	{
		SEQAN_CHECKPOINT;

		//TODO apply diffs of journal to new sequence

		typedef String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
		typedef typename TJournalString::TJournalEntry TJounralEntry;


		JournalTraceDescriptor<TJournalString> traceDescriptor;
		std::stringstream stream;
		stream << journal;
		//can step through trace print and add one journal
		//modify trace descriptor to hold entire insertion information
		appendValue(getTrace(traceDescriptor), TJounralEntry(SOURCE_PATCH, 0, 0, length(journal)));
		append(getInsertionBuffer(traceDescriptor), stream.str());
		applyOperations(journal, host(refJournal), traceDescriptor);
	}
}

#endif /* STRING_SET_JOURNALED_JOIN_SIMPLE_H_ */
