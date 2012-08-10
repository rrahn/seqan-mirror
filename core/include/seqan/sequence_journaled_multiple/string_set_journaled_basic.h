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
  20.12.2010
  ============================================================================
  Implements the basic string set for journaled sequences.
  ==========================================================================*/

#ifndef STRING_SET_JOURNALED_BASIC_H_
#define STRING_SET_JOURNALED_BASIC_H_

namespace seqan {

	struct JournalSetBasic_{};
	typedef Tag<JournalSetBasic_> const JournalSetBasic;

	template <typename TSet>
	struct StringSetGroups
	{
		typedef Nothing Type;
	};


	template <typename TValue>
	class StringSet< TValue, Owner<JournalSetBasic> >
	{
	public:
		typedef typename Position<StringSet>::Type TPosition_;
		typedef String<TValue>	TStrings_;
		typedef typename StringSetLimits<StringSet>::Type TLimits_;
		typedef typename StringSetGroups<StringSet>::Type TGroups_;
		typedef String<TPosition_> TGroupMap_;

		TStrings_	strings;
		TGroups_	groups;
		TGroupMap_  groupMap;

		TLimits_	limits;
		bool		limitsValid;
		bool		groupsValid;

		StringSet() : limitsValid(true), groupsValid(true)
		{
			SEQAN_CHECKPOINT;
			appendValue(limits,0);
		}

		~StringSet()
		{
			SEQAN_CHECKPOINT;
		}

		template <typename TPosition>
		inline typename Reference<StringSet>::Type
		operator[](TPosition const pos)
		{
			return value(*this, pos);
		}

		template <typename TPosition>
		inline typename Reference<StringSet const>::Type
		operator[](TPosition const pos) const
		{
			return value(*this, pos);
		}
	};

	//_____________________Meta________________

	template <typename TValue>
	struct StringSetGroups<StringSet<TValue, Owner<JournalSetBasic> > >
	{
		typedef Group<StringSet<TValue, Owner<JournalSetBasic> >, GroupBasic > TGroup_;
		typedef String<TGroup_> Type;
	};

	template <typename TValue>
	struct StringSetGroups<StringSet<TValue, Owner<JournalSetBasic> > const >
		   : StringSetGroups<StringSet<TValue, Owner<JournalSetBasic> > >{};

	template <typename TValue>
	struct Value<StringSet<TValue, Owner<JournalSetBasic> > >
	{
		typedef TValue Type;
	};

	template <typename TValue>
	struct Value<StringSet<TValue, Owner<JournalSetBasic> > const >
	{
		typedef TValue const Type;
	};

	template <typename TValue>
	struct Reference<StringSet<TValue, Owner<JournalSetBasic> > >
	{
		typedef TValue & Type;
	};

	template <typename TValue>
	struct Reference<StringSet<TValue, Owner<JournalSetBasic> > const >
	{
		typedef TValue const & Type;
	};

	template <typename TValue>
	struct Position<StringSet<TValue, Owner<JournalSetBasic> > >
	{
		typedef size_t Type;
	};

	template <typename TValue>
	struct Position<StringSet<TValue, Owner<JournalSetBasic> > const >
		   : Position<StringSet<TValue, Owner<JournalSetBasic> > >{};

	template <typename TValue>
	struct Size<StringSet<TValue, Owner<JournalSetBasic> > >
	{
		typedef size_t Type;
	};

	template <typename TValue>
	struct Size<StringSet<TValue, Owner<JournalSetBasic> > const >
		   : Size<StringSet<TValue, Owner<JournalSetBasic> > >{};

	//_________________pubic methods_________________

	template <typename TValue, typename TPosition>
	inline typename Reference<StringSet<TValue, Owner<JournalSetBasic> > >::Type
	value(StringSet<TValue, Owner<JournalSetBasic> > & me,
		  TPosition const pos)
	{
		SEQAN_CHECKPOINT;

		me.groupsValid = false;
		return me.strings[pos];
	}

	template <typename TValue, typename TPosition>
	inline typename Reference<StringSet<TValue, Owner<JournalSetBasic> > const >::Type
	value(StringSet<TValue, Owner<JournalSetBasic> > const & me,
		  TPosition const pos)
	{
		SEQAN_CHECKPOINT;

		return me.strings[pos];
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TExpand>
	void appendValue(StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> > & me,
					 String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & newElement,
					 Tag<TExpand> const & tag)
	{
		SEQAN_CHECKPOINT;

		typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
		typedef StringSet<TJournalString, Owner<JournalSetBasic> > TStringSet;
		typedef typename Position<TStringSet>::Type TPosition;
		typedef typename Host<TJournalString>::Type THost;
		typedef typename Value<typename StringSetGroups<TStringSet>::Type >::Type TGroup;

		if(!me.groupsValid)
		{
			resolveGroupMap_(me);
		}
		TGroup group;
		TPosition currPosition = length(me);

		appendValue(me.strings, newElement, tag);

		//create and add group to set

		join(group,currPosition);
		swap(group, currPosition);
		appendValue(me.groups, group, tag);
		appendValue(me.groupMap, currPosition, tag);
	}


	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TExpand>
	void
	appendValue(StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> > & me,
			String<TValue, THostSpec> const & newElement,
			Tag<TExpand> const & tag)
	{
		SEQAN_CHECKPOINT;

		typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
		typedef StringSet<TJournalString, Owner<JournalSetBasic> > TStringSet;
		typedef typename Position<TStringSet>::Type TPosition;
		typedef typename Host<TJournalString>::Type THost;
		typedef typename Value<typename StringSetGroups<TStringSet>::Type >::Type TGroup;

		if (!me.groupsValid)
		{
			resolveGroupMap_(me);
		}

		TJournalString jrn;
		setHost(jrn, newElement);
		TGroup group;

		join(group, length(me));
		swap(group, length(me));
		appendValue(me.strings, jrn, tag);
		appendValue(me.groupMap, length(me.groups), tag);
		appendValue(me.groups, group, tag);
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TExpand>
	void
	assignValue(StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> > & me,
				String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & newValue,
				typename Position<StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> > >::Type const & pos)
	{
		SEQAN_CHECKPOINT;

		//check for valid string position
		SEQAN_ASSERT_LEQ(pos, length(me));

		typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
		typedef StringSet<TJournalString, Owner<JournalSetBasic> > TStringSet;
		typedef typename Position<TStringSet>::Type TPosition;
		typedef typename Host<TJournalString>::Type THost;
		typedef typename Value<typename StringSetGroups<TStringSet>::Type >::Type TGroup;

		if (!me.groupsValid)
		{
			resolveGroupMap_(me);
		}


		value(me, pos) == newValue;

		if(&host(me[getGroupReference(me.groups[me.groupMap[pos]])]) != &host(newValue))
		{
			TPosition relatedGroupPos = findGroup_(me, pos);
			if (relatedGroupPos == maxValue<TPosition>())
			{//create new group and append it
				TGroup newGroup;

				join(newGroup, pos);
				swap(newGroup, pos);
				relatedGroupPos = length(me.groups);
				appendValue(me.groups, newGroup, Generous());
			}
			//disjoin value from old group
			disjoin(me.groups[me.groupMap[pos]], pos);
			//join to new group
			join(me.groups[relatedGroupPos], pos);
			//update groupMap
			me.groupMap[pos] == relatedGroupPos;
		}
		 //TODO rmaerker: assert consistency check
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TExpand>
	void
	assignValue(StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> > & me,
				String<TValue, THostSpec> & newValue,
				typename Position<StringSet<String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> > >::Type const & pos)
	{
		SEQAN_CHECKPOINT;
		//check for valid string position
		SEQAN_ASSERT_LEQ(pos, length(me));

		typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
		typedef StringSet<TJournalString, Owner<JournalSetBasic> > TStringSet;
		typedef typename Position<TStringSet>::Type TPosition;
		typedef typename Host<TJournalString>::Type THost;
		typedef typename Value<typename StringSetGroups<TStringSet>::Type >::Type TGroup;

		if (!me.groupsValid)
		{
			resolveGroupMap_(me);
		}

		//first create journal of new assigned sequence
		TJournalString jrn;
		setHost(jrn, newValue);

		value(me, pos) == jrn;
		//TODO rmaerker: what happens with the group if only one sequence is contained
		//TODO rmaerker: what happens with the group ref if it is reassigened with another value? - join new assigned sequence to the group
		//can not have a group that refers to a jrn that does not exist anymore within the set
		if(&host(me[getGroupReference(me.groups[me.groupMap[pos]])]) != &host(newValue))
		{
			TPosition relatedGroupPos = findGroup_(me, pos);
			if (relatedGroupPos == maxValue<TPosition>())
			{//create new group and append it
				TGroup newGroup;

				join(newGroup, pos);
				swap(newGroup, pos);
				relatedGroupPos = length(me.groups);
				appendValue(me.groups, newGroup, Generous());
			}
			//disjoin value from old group
			disjoin(me.groups[me.groupMap[pos]], pos);
			//join to new group
			join(me.groups[relatedGroupPos], pos);
			//update groupMap
			me.groupMap[pos] == relatedGroupPos;
		}
		 //TODO rmaerker: assert resolved goupMap
	}

	template <typename TValue, typename TSize>
	inline void
	resize(StringSet<TValue, Owner<JournalSetBasic> > & me,
			TSize const & newSize)
	{
		SEQAN_CHECKPOINT;

		resize(me.strings, newSize);
		resize(me.groupMap, newSize);
		resize(me.groups, newSize);
	}

	template <typename TValue>
	inline void
	clear(StringSet<TValue, Owner<JournalSetBasic> > & me)
	{
		SEQAN_CHECKPOINT;

		clear(me.strings);
		clear(me.groups);
		clear(me.groupMap);
		clear(me.limits);
		me.limitsValid = true;
		me.groupsValid = true;
	}

	//___________________internal methods____________________

	template <typename TValue, typename TSetSpec>
	inline typename Position<StringSet<TValue, Owner<TSetSpec> > const >::Type
	findGroup_(StringSet<TValue, Owner<TSetSpec> > const & me,
				typename Position<StringSet<TValue, Owner<TSetSpec> > const >::Type const & pos)
	{
		SEQAN_CHECKPOINT;

		typedef typename Position<StringSet<TValue, Owner<TSetSpec> > const >::Type TPosition;

		TValue const &jrnAdr = value(me, pos);

		for (TPosition i = 0u; i < length(getGroups(me)); ++i)
		{
			TValue const &refAdr = me[getGroupReference(me.groups[i])];
			if (&host(jrnAdr) == &host(refAdr)) {
				return i;
			}
		}
		return maxValue<TPosition>();
	}

	template <typename TValue>
	inline void
	resolveGroupMap_(StringSet<TValue, Owner<JournalSetBasic> > & me)
	{
		SEQAN_CHECKPOINT;
		typedef typename Host<TValue>::Type THost;
		typedef StringSet<TValue, Owner<JournalSetBasic> > TStringSet;
		typedef typename Position<TStringSet>::Type TPosition;
		typedef typename Value<typename StringSetGroups<TStringSet>::Type >::Type TGroup;

		for (TPosition i = 0u; i < length(me.groupMap);++i)
		{
			TPosition groupPos = value(me.groupMap, i);
			TPosition groupRefPos = getGroupReference(value(me.groups, groupPos));
			THost const &jrnHostPtr = host(value(me.strings,i));
			THost const &grpHostPtr = host(value(me.strings, groupRefPos));
			//journal doesn't belong to group used in groupMap.
			//changing the journal in the set by assigning a new journal
			if(jrnHostPtr != grpHostPtr)
			{
				//maybe journal is assigned to other group - find group
				TPosition relatedGroupPos = findGroup_(me, i);
				//no group exist for journal - add new one
				if (relatedGroupPos == maxValue<TPosition>())
				{
					TGroup newGroup;

					join(newGroup, i);
					swap(newGroup, i);
					relatedGroupPos = length(me.groups);
					appendValue(me.groups, newGroup, Generous());
				}
				//disjoin value from old group
				disjoin(me.groups[me.groupMap[i]], i);
				//join to new group
				join(me.groups[relatedGroupPos], i);
				//update groupMap
				me.groupMap[i] = relatedGroupPos;
			}
		}
		me.groupsValid = true;
	}

	//_______________________ GETTER _________________________________

	template <typename TValue>
	inline String<TValue> &
	getStrings(StringSet<TValue, Owner<JournalSetBasic> > & me)
	{
		SEQAN_CHECKPOINT;
		return me.strings;
	}

	template <typename TValue>
	inline String<TValue> const &
	getStrings(StringSet<TValue, Owner<JournalSetBasic> > const & me)
	{
		SEQAN_CHECKPOINT;
		return me.strings;
	}

	template <typename TValue>
	inline typename StringSetGroups<StringSet<TValue, Owner<JournalSetBasic> > >::Type &
	getGroups(StringSet<TValue, Owner<JournalSetBasic> > & me)
	{
		SEQAN_CHECKPOINT;
		return me.groups;
	}

	template <typename TValue>
	inline typename StringSetGroups<StringSet<TValue, Owner<JournalSetBasic> > >::Type const &
	getGroups(StringSet<TValue, Owner<JournalSetBasic> > const & me)
	{
		SEQAN_CHECKPOINT;
		return me.groups;
	}

	template <typename TValue>
	inline String<typename Position<StringSet<TValue, Owner<JournalSetBasic> > >::Type> &
	getGroupMap(StringSet<TValue, Owner<JournalSetBasic> > & me)
	{
		SEQAN_CHECKPOINT;
		return me.groupMap;
	}

	template <typename TValue>
	inline String<typename Position<StringSet<TValue, Owner<JournalSetBasic> > >::Type> const &
	getGroupMap(StringSet<TValue, Owner<JournalSetBasic> > const & me)
	{
		SEQAN_CHECKPOINT;
		return me.groupMap;
	}


//	template <typename TValue>
//	class StringSet< TValue, Owner<JournalSetBasic> > {
//	public:
//		//______typedefs_________________________________________________________________________
//
//		typedef String< TValue >  TJournals_;
//		typedef JournalGroup<StringSet> TGroup_;
//		typedef String<TGroup_>		    TGroupString_;
//		typedef typename StringSetLimits<StringSet>::Type	TLimits;
//		typedef typename Concatenator<StringSet>::Type		TConcatenator;
//
//		//______member variables_________________________________________________________________
//
//		TJournals_			strings;		//contains all assigned journals
//		TGroupString_		groups;			//contains the groups the journals are assigned to
//		String<int>		   groupMap; 		//lists the group each journal is assigned to.
//
//		//internal variables
//		TLimits			limits;				// needed for global string set class ... see sequence/multiple_sequence.h
//		bool			limitsValid;		// is true if limits contains the cumulative sum of the sequence lengths
//		TConcatenator	concat;
//
//		//_____constructor_______________________________________________________________________
//
//		//default ctor
//		StringSet(): limitsValid(true) {
//			SEQAN_CHECKPOINT
//				appendValue(limits, 0);
//		}
//
//		//copy-ctor
//		StringSet(StringSet < TValue, Owner< JournalSetBasic > > const & other) {
//			SEQAN_CHECKPOINT;
//			if (*this != other) {
//				//resize(*this, length(other));
//				//for (unsigned int i = 0; i < length(other.journals); ++i)
//				//{
//				//	journals[i] = other.journals[i];  //copies the journals in the new set
//				//}
//				//for (unsigned int i = 0; i < length(other.groups); ++i)
//				//{
//				//	groups[i] = other.groups[i];  //copies the groups in the new set
//				//}
//				limits = other.limits;
//			}
//		}
//
//		~StringSet() {
//			SEQAN_CHECKPOINT;
//				//nothing to do
//		}
//
//		//_____operator_______________________________________________________________________
//
//		template <typename TPos>
//		inline typename Reference<StringSet>::Type
//			operator [] (TPos pos)
//		{
//			SEQAN_CHECKPOINT;
//			return value(*this, pos);
//		}
//
//		template <typename TPos>
//		inline typename Reference<StringSet const>::Type
//			operator [] (TPos pos) const
//		{
//			SEQAN_CHECKPOINT;
//			return value(*this, pos);
//		}
//	};



}

#endif /* STRING_SET_JOURNALED_BASIC_H_ */
