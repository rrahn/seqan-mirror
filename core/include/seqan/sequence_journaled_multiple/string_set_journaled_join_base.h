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
  Implements the basic interface to choose the correct join methods and to
  generate journal trees out of the alignment results
  ==========================================================================*/

#ifndef STRING_SET_JOURNALED_JOIN_BASE_H_
#define STRING_SET_JOURNALED_JOIN_BASE_H_

namespace SEQAN_NAMESPACE_MAIN
{

	//___________________Meta_____________________________
	struct SimpleJoin_;
	typedef Tag<SimpleJoin_> const SimpleJoin;

	struct MaxCompression_;
	typedef Tag<MaxCompression_> const MaxCompression;

	/**
	 * Returns the default score for the small join method.
	 */
	template <typename TValue>
	inline Score<TValue, ExtendedScore>
	defaultMaxCompressionScore()
	{
		SEQAN_CHECKPOINT;

		return Score<TValue, ExtendedScore>(0,-(maxValue<TValue>() >> 1), -25,-1,-12,0);
	}


	//TODO rmaerker: structure of the join operations
	//TODO rmaerker: return achieved compression ratio
	template <typename TValue, typename TGroupId, typename TJournalId, typename TScoreValue, typename TScoreSpec, typename TStrategy>
	void
	join(StringSet<TValue, Owner<JournalSetBasic> > & me,
			  TGroupId groupId,
			  TJournalId journalId,
			  Score<TScoreValue, TScoreSpec> const & score = defaultMaxCompressionScore<TScoreValue>(),
			  TStrategy const & strategy = MaxCompression())
	{
		SEQAN_CHECKPOINT;

		SEQAN_ASSERT_LEQ(groupId, length(me));
		SEQAN_ASSERT_LEQ(journalId, length(me));

		typedef StringSet<TValue, Dependent<Tight> > TAlignSet;

		//rmaerker distinguish join operation
		if (!appliesTo(me[journalId]), host(me[getGroupReference(getGroups(me)[groupId])]))
		{
			internalJoin_( me[journalId], me[getGroupReference(getGroups(me)[groupId])], score, strategy);
		}

		updateGroup(me, groupId, journalId);
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	inline bool
	appliesTo(String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & journal,
			  String <TValue, THostSpec> const & reference)
	{
		SEQAN_CHECKPOINT;
		return (&host(journal) == &reference);
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	inline void
	applyOperations(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journal,
			     String<TValue, THostSpec> const & newHost,
			     JournalTraceDescriptor<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > > const & traceDescr)
	{
		SEQAN_CHECKPOINT;

		setValue(journal._holder, newHost);
		reinit(journal._journalEntries, length(newHost));
		constructAndSetJournalTree(journal, traceDescr);
	}

	/**
	 * Updates the group information for the groupId and the journalId within the
	 * given set, such that the journalId is disjoined from its old group and
	 * joined to the given group.
	 */
	template <typename TValue, typename TGroupId, typename TJournalId>
	inline void
	updateGroup(StringSet<TValue, Owner<JournalSetBasic> > & me,
			      TGroupId groupId,
			      TJournalId journalId)
	{
		SEQAN_CHECKPOINT;

		SEQAN_ASSERT_LEQ(groupId, length(me));
		SEQAN_ASSERT_LEQ(journalId, length(me));

		disjoin(me.groups[me.groupMap[journalId]], journalId);
		join(groupId,journalId);
		me.groupMap[journalId] = groupId;
		me.groupsValid = true;
		//TODO rmaerker: implement check consistent group map
//		SEQAN_ASSERT_TRUE(checkGroupMap_(me));
	}

}
#endif /* STRING_SET_JOURNALED_JOIN_BASE_H_ */
