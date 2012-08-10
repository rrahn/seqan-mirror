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
  Implements the interface for the join operation which achieves the highest
  compression ratio between two sequences
  ==========================================================================*/

#ifndef STRING_SET_JOURNALED_JOIN_MAX_COMPRESSION_H_
#define STRING_SET_JOURNALED_JOIN_MAX_COMPRESSION_H_

namespace SEQAN_NAMESPACE_MAIN
{
	/**
	 * Computes the highest compression between two given sequences.
	 *
	 * Returns the compression ratio obtained by the algorithm.
	 */
	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TScore>
	float
	internalJoin_(String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journal,
				  String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & refJournal,
				  TScore const &  score,
				  MaxCompression const &)
	{
		SEQAN_CHECKPOINT;

		typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
		typedef typename Position<TJournalString>::Type TPosition;
		typedef typename Host<TJournalString>::Type THost;
		typedef typename Size<TJournalString>::Type TSize;
		typedef typename JournalType<TJournalString>::Type TJournalTree;
		typedef typename Value<TScore>::Type TScoreValue;

		typedef StringSet<TJournalString, Dependent<> > TAlignSet;
		TSize oldLength = length(journal);

		TJournalString refOnlyJournal;
		setHost(refOnlyJournal, host(refJournal));

		TAlignSet alignSet;
		appendValue(alignSet, journal);
		appendValue(alignSet, refOnlyJournal);

		//assure we compare only to the ref, not to the journal of the ref, which might has a tre
		String<unsigned char> trace;
		unsigned char initialDir;
		TScoreValue overallMaxValue[2];
		TSize overallMaxIndex[2];

		JournalTraceDescriptor<TJournalString> traceDescriptor;

		_alignGotoh(trace,alignSet, score, overallMaxValue, overallMaxIndex, initialDir, AlignConfig<>());
		_alignGotohTrace(traceDescriptor, alignSet, trace, overallMaxIndex, initialDir);

		std::cout << traceDescriptor << std::endl;
		applyOperations(journal, host(refJournal), traceDescriptor);
		return oldLength/length(traceDescriptor.insertionBuffer_);
	}

}

#endif /* STRING_SET_JOURNALED_JOIN_MAX_COMPRESSION_H_ */
