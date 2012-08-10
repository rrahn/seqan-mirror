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
 Implements the group class used to group journals
 ==========================================================================*/

#ifndef STRING_SET_JOURNALED_GROUP_H_
#define STRING_SET_JOURNALED_GROUP_H_

namespace seqan
{

	struct GroupBasic_;
	typedef Tag<GroupBasic_> const GroupBasic;

	template<typename TSet, typename TGroupSpec>
	class Group
	{

	};

	template<typename TSet>
	class Group<TSet, GroupBasic>
	{
		public:

			typedef typename Position<TSet>::Type 	TPosition;
			typedef String<TPosition> 				TMembers;
			typedef String<char> 					TLabel;

			//_______________________

			TPosition reference_;
			TMembers members_;
			TLabel label_;

			Group() :
				reference_(MaxValue<TPosition>::VALUE)
			{
				SEQAN_CHECKPOINT;
			}

			Group(Group const & other) :
				reference_(other.reference_), members_(other.members_), label_(
						other.label_)
			{
				SEQAN_CHECKPOINT;
			}

			~Group()
			{
				SEQAN_CHECKPOINT;
			}

			//________________

			template<typename TPosition>
			inline typename Reference<Group>::Type operator[](
					TPosition const pos)
			{
				SEQAN_CHECKPOINT;
				return value(*this, pos);
			}

			template<typename TPosition>
			inline typename Reference<Group const>::Type operator[](
					TPosition const pos) const
			{
				SEQAN_CHECKPOINT;
				return value(*this, pos);
			}

			Group &
			operator=(Group const & other)
			{
				SEQAN_CHECKPOINT;
				if (*this != other)
				{
					setGroupLabel(*this, other.label_);
					setGroupReference(*this, other.reference_);
					setGroupMembers(*this, other.members_);
				}
				return *this;
			}
	};

	//____________________ META _________________________

	template<typename TSet, typename TGroupSpec>
	struct Value<Group<TSet, TGroupSpec> >
	{
			typedef typename Position<TSet>::Type Type;
	};

	template<typename TSet, typename TGroupSpec>
	struct Value<Group<TSet, TGroupSpec> const>
	{
			typedef typename Position<TSet>::Type const Type;
	};

	template<typename TSet, typename TGroupSpec>
	struct Reference<Group<TSet, TGroupSpec> >
	{
			typedef typename Value<Group<TSet, GroupBasic> >::Type & Type;
	};

	template<typename TSet, typename TGroupSpec>
	struct Reference<Group<TSet, TGroupSpec> const>
	{
			typedef typename Value<Group<TSet, GroupBasic> >::Type const & Type;
	};

	template < typename TSet, typename TGroupSpec>
	struct Size<Group<TSet, TGroupSpec> > {
		typedef typename Size< TSet >::Type Type;
	};

	template < typename TSet, typename TGroupSpec>
	struct Size<Group<TSet, TGroupSpec> const>
		: Size<Group<TSet, TGroupSpec> >{};

	template < typename TSet, typename TGroupSpec>
	struct Position<Group<TSet, TGroupSpec> > {
		typedef typename Size< TSet >::Type Type;
	};

	template < typename TSet, typename TGroupSpec>
	struct Position<Group<TSet, TGroupSpec> const>
		: Position<Group<TSet, TGroupSpec> >{};

	template < typename TGroup>
	struct GroupMembers{};

	template < typename TSet, typename TGroupSpec >
	struct GroupMembers< Group<TSet, TGroupSpec> >
	{
		typedef typename Value < Group < TSet, TGroupSpec> >::Type TValue_;
		typedef String < TValue_ > Type;
	};

	template < typename TSet, typename TGroupSpec >
	struct GroupMembers< Group<TSet, TGroupSpec> const >
		//: JournalGroupMembers< JournalGroup < TJournalSet, JournalGroupSpec < TGroupSpec > > >
	{
		typedef typename Value < Group < TSet, TGroupSpec> >::Type TValue_;
		typedef String < TValue_ > const Type;
	};


	template < typename TGroup >
	struct GroupLabel{};

	template < typename TSet, typename TGroupSpec >
	struct GroupLabel< Group < TSet, TGroupSpec> >
	{
		typedef String < char > Type;
	};

	template < typename TSet, typename TGroupSpec >
	struct GroupLabel< Group < TSet, TGroupSpec > const >
	{
		typedef String < char > const Type;
	};

	//____________________ Function _________________________

	template<typename TSet, typename TGroupSpec, typename TPosition>
	inline typename Reference<Group<TSet, TGroupSpec> >::Type
	value(Group<TSet,TGroupSpec> & me,
		  TPosition const & pos)
	{
		SEQAN_CHECKPOINT;
		return me.members_[pos];
	}

	template<typename TSet, typename TGroupSpec, typename TPosition>
	inline typename Reference<Group<TSet, TGroupSpec> const>::Type
	value(Group<TSet, TGroupSpec> const & me,
		  TPosition const & pos)
	{
		SEQAN_CHECKPOINT;
		return me.members_[pos];
	}

	/*
	 *	Sets the label to the group
	 */
	template<typename TJournalSet, typename TGroupSpec>
	inline void
	setGroupLabel(Group<TJournalSet, TGroupSpec> & group,
				  String<char> const & newLabel)
	{
		SEQAN_CHECKPOINT;
		group.label_ = newLabel;
	}

	/*
	 *	Returns the label of the group
	 */
	template<typename TJournalSet, typename TGroupSpec>
	inline String<char> &
	getGroupLabel(Group<TJournalSet, TGroupSpec> & group)
	{
		SEQAN_CHECKPOINT;
		return group.label_;
	}

	/*
	 *	Returns the label of the group
	 */
	template<typename TJournalSet, typename TGroupSpec>
	inline String<char> const &
	getGroupLabel(Group<TJournalSet, TGroupSpec> const & group)
	{
		SEQAN_CHECKPOINT;
		return group.label_;
	}

	template<typename TJournalSet, typename TGroupSpec>
	inline void
	setGroupReference(Group<TJournalSet, TGroupSpec> & group,
					  typename Position<Group<TJournalSet, TGroupSpec> >::Type const & refId)
	{
		SEQAN_CHECKPOINT;
		group.reference_ = refId;
	}

	template<typename TJournalSet, typename TGroupSpec>
	inline typename Position<Group<TJournalSet, TGroupSpec> >::Type getGroupReference(
			Group<TJournalSet, TGroupSpec> & group)
	{
		SEQAN_CHECKPOINT;
		return group.reference_;
	}

	template<typename TJournalSet, typename TGroupSpec>
	inline typename Position<Group<TJournalSet, TGroupSpec> const>::Type getGroupReference(
			Group<TJournalSet, TGroupSpec> const & group)
	{
		SEQAN_CHECKPOINT;
		return group.reference_;
	}

	template<typename TJournalSet, typename TGroupSpec>
	inline void setGroupMembers(Group<TJournalSet, TGroupSpec > & group,
			String<typename Position<Group<TJournalSet, TGroupSpec> >::Type> newMembers)
	{
		SEQAN_CHECKPOINT;
		resize(group.members_, length(newMembers));
		for (unsigned int i = 0; i < length(newMembers); ++i)
		{
			group.members_[i] = newMembers[i];
		}
	}

	template<typename TJournalSet, typename TGroupSpec>
	inline typename GroupMembers<Group<TJournalSet,TGroupSpec> >::Type &
	getGroupMembers(Group<TJournalSet, TGroupSpec> & group)
	{
		SEQAN_CHECKPOINT;
		return group.members_;
	}

	template<typename TJournalSet, typename TGroupSpec>
	inline typename GroupMembers<Group<TJournalSet,TGroupSpec> const>::Type &
	getGroupMembers(Group<TJournalSet, TGroupSpec> const & group)
	{
		SEQAN_CHECKPOINT;
		return group.members_;
	}

	/*
	 * Checks two journal groups for equality.
	 */
	template<typename TJournalSet, typename TGroupSpec>
	inline bool operator==(Group<TJournalSet, TGroupSpec> const & lObject,
			Group<TJournalSet, TGroupSpec> const & rObject)
	{
		SEQAN_CHECKPOINT;

		if (lObject.label_ == rObject.label_)
		{
			if (lObject.reference_ == rObject.reference_)
			{
				if (lObject.members_ == rObject.members_)
				{
					return true;
				}
			}
		}
		return false;
	}

	/*
	 * Checks two journal groups for inequality.
	 */
	template<typename TJournalSet, typename TGroupSpec>
	inline bool
	operator!=(Group<TJournalSet, TGroupSpec> const & lObject,
			Group<TJournalSet, TGroupSpec> const & rObject)
	{
		SEQAN_CHECKPOINT;
		return !(lObject == rObject);
	}

	/**
	 *Clears the container.
	 */
	template < typename TSet, typename TGroupSpec >
	inline void
	clear(Group < TSet, TGroupSpec > & me)
	{
		SEQAN_CHECKPOINT
		typedef typename Position<Group<TSet, TGroupSpec> >::Type TPosition;

		setGroupReference(me,MaxValue<TPosition>::VALUE);
		setGroupLabel(me, "");
		resize(getGroupMembers(me), 0);
	}

	template < typename TSet, typename TGroupSpec, typename TExpand >
	inline void
	appendValue(Group <TSet, TGroupSpec> & me,
				typename Position <Group<TSet, TGroupSpec> >::Type const & memberId,
				Tag < TExpand > const & tag)
	{
		SEQAN_CHECKPOINT;
		appendValue(getGroupMembers(me), memberId, tag);
	}

	/*
	 *	Returns true if the given journal position is equal to the reference id of the group.
	 */
	template < typename TSet, typename TGroupSpec, typename TPosition >
	inline bool
	isReference(Group < TSet, TGroupSpec> & me,
					TPosition const & pos)
	{
		SEQAN_CHECKPOINT;

		if (pos == getGroupReference(me))
		{
			return true;
		}
		return false;
	}

	template <typename TJournalGroup, typename TPosition >
	inline typename Size < TJournalGroup >::Type
	findPosition_(TJournalGroup const & me,
				 TPosition const & pos)
	{
		SEQAN_CHECKPOINT;
		for (TPosition i = 0; i < (TPosition) length(me); ++i)
		{
			if (pos == value(getGroupMembers(me),i)){
				return i;
			}
		}
		return length(me);
	}

	/*
	 *	Returns the length of the group
	 */
	template < typename TSet, typename TGroupSpec >
	inline typename Size < Group < TSet, TGroupSpec > const >::Type
	length( Group < TSet, TGroupSpec > const & me)
	{
		SEQAN_CHECKPOINT;
		return length(me.members_);
	}

	/*
	 *	Changes the id that represents the reference of the group to the given id.
	 */
	template < typename TSet, typename TGroupSpec >
	inline bool
	swap( Group < TSet, TGroupSpec > & me,
		 typename Position<Group<TSet, TGroupSpec> >::Type const & pos)
	{
		SEQAN_CHECKPOINT;
		SEQAN_CHECKPOINT;
		if (isMember(me, pos))
		{
			setGroupReference(me, pos);
			return true;
		}
		return false;
	}

	/*
	 *	Renames the group.
	 */
	template < typename TSet, typename TGroupSpec >
	inline void
	rename(Group < TSet, TGroupSpec > & me,
				String <char> const & newLabel)
	{
		SEQAN_CHECKPOINT;
		setGroupLabel(me, newLabel);
	}

	template < typename TJournalGroup, typename TPosition >
	inline bool
	isMember(TJournalGroup const & me,
			TPosition const & pos)
	{
		SEQAN_CHECKPOINT;
		for (TPosition i = 0; i < (TPosition) length(me);++i)
		{
			if (pos == value(getGroupMembers(me), i)) {
				return true;
			}
		}
		return false;
	}

	template < typename TJournalGroup, typename TPosition >
	inline bool
	join(TJournalGroup & me,
		TPosition const & pos)
	{
		SEQAN_CHECKPOINT;
		return join(me, pos, Generous());
	}

	template < typename TJournalGroup, typename TPosition, typename TExpand >
	inline bool
	join(TJournalGroup & me,
		TPosition const & pos,
		Tag<TExpand> const & tag)
	{
		SEQAN_CHECKPOINT;
		if (isMember(me, pos)) {
			return false;
		}
		appendValue(me, pos, tag);
		return true;
	}

	template <typename TSet, typename TGroupSpec, typename TPosition>
	inline void
	erase(Group < TSet, TGroupSpec > & me,
			 TPosition const & pos)
	{
		SEQAN_CHECKPOINT;
		erase(getGroupMembers(me), pos);
	}

	//TODO rmaerker: either allow swaping the reference or not
	template < typename TGroup, typename TPosition >
	inline bool
	disjoin(TGroup & me,
			TPosition const & pos)
	{
		SEQAN_CHECKPOINT;

		if (isReference(me, pos) || !isMember(me, pos))
		{
			return false;
		}
		TPosition delPos = findPosition_(me, pos);
		SEQAN_ASSERT_TRUE(delPos < length(me));
		erase(me, delPos);
		return true;
	}
}

#endif /* STRING_SET_JOURNALED_GROUP_H_ */
