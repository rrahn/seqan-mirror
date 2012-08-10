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
  06.01.2011
  ============================================================================
  Defines the utility methods for the string set
  ==========================================================================*/

#ifndef STRING_SET_JOURNALED_UTIL_H_
#define STRING_SET_JOURNALED_UTIL_H_

namespace seqan{



	enum NodeType
	{
		NODE_NULL,
		NODE_INSERTION,
		NODE_DELETION
	};

	template <typename TPosition_, typename TSize_, typename TInsertBuff_>
	struct InfoNode
	{
			typedef TPosition_ TPosition;
			typedef TSize_ TSize;
			typedef TInsertBuff_ TInsertBuff;

			NodeType	nodeType_;
			TPosition	begin_;
			TSize 		length_;
			TInsertBuff insertionBuff_;

			InfoNode() :
				nodeType_(NODE_NULL), begin_(0), length_(0)
			{
				SEQAN_CHECKPOINT;
			}

			InfoNode(InfoNode const & other) :
				nodeType_(other.nodeType_), begin_(other.begin_), length_(other.length_), insertionBuff_(other.insertionBuff_)
			{
				SEQAN_CHECKPOINT;
			}

			~InfoNode()
			{
				SEQAN_CHECKPOINT;
			}
	};


	//___________________public methods___________________

	template <typename TValue, typename TPosition, typename TScore, typename TStrategy>
	void
	merge(StringSet<TValue, Owner<JournalSetBasic> > & me,
		TPosition const & targetGrpPos,
		TPosition const & sourceGrpPos,
		TScore const & score = Score<TPosition>(),
		TStrategy const & tag = MaxCompression())
	{
		SEQAN_CHECKPOINT;

		typedef StringSet<TValue, Owner<JournalSetBasic> > TJournalSet;
		typedef typename Value< typename StringSetGroups<TJournalSet>::Type >::Type TGroup;

		TGroup target = me.groups[targetGrpPos];
		TGroup source = me.groups[sourceGrpPos];

		if (!appliesTo(me[getGroupReference(target)], me[getGroupReference(source)]))
		{
			internalMerge_(me, target, source, score, tag);
		}

		//completely remove the group
	}

	template <typename TValue, typename TPosition, typename TScore, typename TStrategy>
	void
	join(StringSet<TValue, Owner<JournalSetBasic> > & me,
		 TPosition & groupPos,
		 TPosition & jrnPos,
		 TScore const & score = Score<TPosition>(),
		 TStrategy const & tag = MaxCompression())
	{
		SEQAN_CHECKPOINT;

		typedef StringSet<TValue, Owner<JournalSetBasic> > TJournalSet;
		typedef typename Value< typename StringSetGroups<TJournalSet>::Type >::Type TGroup;

		TGroup targetGroup = me.groups[groupPos];
		TValue target = me[getGroupReference(targetGroup)];
		TValue source = me[jrnPos];

		if (!appliesTo(target, source))
		{
			internalJoin_(target, source, score, tag);
		}
		//chek because u could clear a group... if only the reference is in
		//what happens if the reference should be disjoined?
		disjoin(me.groups[me.groupMap[jrnPos]], jrnPos);
		join(targetGroup,jrnPos);
		me.groupMap[jrnPos] = groupPos;
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	inline bool
	appliesTo(String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & target,
			String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & source)
	{
		SEQAN_CHECKPOINT;
		return (&host(target) == &host(source));
	}


	//___________________internal methods___________________

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec>
	void
	adaptGroupMembers_(StringSet< String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> > & set,
				     String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & source)
	{
		SEQAN_CHECKPOINT;

		typedef String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
		typedef StringSet< TJournalString, Owner<JournalSetBasic> > TJournalSet;
		typedef typename Value< typename StringSetGroups<TJournalSet>::Type >::Type TGroup;
		typedef typename Position<TGroup>::Type TPosition;

		if (!isFlat(set[getGroupReference(source)]))
		{
			for (TPosition i = 0; i < length(source); ++i)
			{
				if (getGroupMembers(source)[i] != getGroupReference(source))
				{
					//TODO rmaerker: apply new changes to the journal
					//find a way to obtain
					// how can i do this best
					//step through journal and
				}
			}
		}
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TScore, typename TStrategy>
	void
	internalMerge_(StringSet< String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> > & set,
				   String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & target,
				   String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & source,
				   TScore const & score,
				   TStrategy const & tag)
	{
		SEQAN_CHECKPOINT;
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TScore>
	void
	internalMerge_(StringSet< String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> > & set,
				   Group <StringSet< String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> >, GroupBasic >  const & target,
				   Group <StringSet< String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> >, GroupBasic > & source,
				   TScore const & /*score*/,
				   SimpleJoin const & /*strategy*/)
	{
		SEQAN_CHECKPOINT;

		typedef String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
		typedef StringSet<TJournalString, Owner<JournalSetBasic> > TJournalSet;
		typedef typename Value< typename StringSetGroups<TJournalSet>::Type >::Type TGroup;
		typedef typename Position<TGroup>::Type TGrpPos;
		typedef typename Host<TJournalString>::Type THost;
		typedef typename JournalType<TJournalString>::Type TJournalType;
		typedef typename Iterator<TJournalType>::Type TIterator;

		TJournalString targetRef = set[getGroupReference(target)];

		for (TGrpPos i = 0; i < length(source); ++i)
		{
			internalJoin_(targetRef, set[source[i]], TScore(), SimpleJoin());
		}
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TScore>
	void
	internalMerge_(StringSet< String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> > & set,
			Group <StringSet< String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> >, GroupBasic >  const & target,
			Group <StringSet< String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> >, Owner<JournalSetBasic> >, GroupBasic > & source,
			TScore const & /*score*/,
			MaxCompression const & /*strategy*/)
	{
		SEQAN_CHECKPOINT;

		//TODO rmaerker: small join between two reference sequences and than apply changes to all other journals in the source group
		typedef String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
		typedef StringSet<TJournalString, Owner<JournalSetBasic> > TJournalSet;

		typedef typename Value< typename StringSetGroups<TJournalSet>::Type >::Type TGroup;
		typedef typename Position<TGroup>::Type TGrpPos;
		typedef typename Host<TJournalString>::Type THost;
		typedef typename JournalType<TJournalString>::Type TJournalType;
		typedef typename Iterator<TJournalType>::Type TIterator;

		TJournalString targetRef = set[getGroupReference(target)];
		TJournalString sourceRef = set[getGroupReference(source)];

		internalJoin_(targetRef, sourceRef, TScore(), MaxCompression());

		for (TGrpPos i = 0u; i < length(source); ++i)
		{
			if (source[i] != getGroupReference(source))
			{
				synchronizeJournalTrees_(sourceRef, set[source[i]]);
			}
		}


	}

	template <typename TValue, typename THostSpec, typename TBuffSpec>
	void
	synchronizeJournalTrees_(String<TValue, Journaled<THostSpec, SortedArray, TBuffSpec> > & target,
						  String<TValue, Journaled<THostSpec, SortedArray, TBuffSpec> > & source)
	{
		SEQAN_CHECKPOINT;

		typedef String<TValue, Journaled<THostSpec, SortedArray, TBuffSpec> > TJournalString;
		typedef typename Host<TJournalString>::Type THost;
		typedef String<TValue, TBuffSpec> TInsertionBuffer;

		typedef typename JournalType<TJournalString>::Type TJournalEntries;
		typedef typename Value<TJournalEntries>::Type TJournalEntry;
		typedef typename Position<TJournalEntry>::Type TPos;
		typedef typename Size<TJournalEntry>::Type TSize;
		typedef JournalEntryLtByVirtualPos<TPos, TSize> TCmp;

		typedef typename Iterator<TJournalEntries, Standard>::Type TIterator;

		TIterator targetBegin = begin(target._journalEntries);
		TIterator targetEnd = end(target._journalEntries);
		TIterator sourceBegin = begin(source._journalEntries);
		TIterator sourceEnd = end(source._journalEntries);
		TIterator targetIt = targetBegin;

		TJournalEntries * sourceTreePtr = &source._journalEntries;
		//manages the current position within the journal entries of the source
		TJournalEntry currentNode = *targetBegin;

		int shiftPhysicalPos = 0;
		int shiftVirtualPos = 0;
		TPos lastPhysPos = 0;

		//set host information from target to source sequence
		(*sourceTreePtr)._originalStringLength = target._journalEntries._originalStringLength;
		setValue(source._holder, host(target));

		while(targetIt != targetEnd)
		{
			TIterator sourceIt = std::upper_bound(sourceBegin, sourceEnd, currentNode, TCmp());
			//MUST NOT find begin
			SEQAN_ASSERT_TRUE(sourceIt != sourceBegin);
			--sourceIt;

			//current node points behind last node in the source
			if (currentNode.virtualPosition < sourceIt->virtualPosition + sourceIt->length)
			{
				TPos pos = sourceIt - sourceBegin;
				switch(currentNode.segmentSource + (sourceIt->segmentSource << 1))
				{
				case 3u:	//tmp->SOURCE_ORIGINAL & jrn->SOURCE_ORIGINAL
				{
					SEQAN_ASSERT_EQ_MSG(currentNode.virtualPosition, sourceIt->virtualPosition,"Source virtual position <%d> does not match target virtual position <%d>!", currentNode.virtualPosition, sourceIt->virtualPosition);
					//update physical position first
					sourceIt->physicalPosition += shiftPhysicalPos;

					//case 1: deleted part of target node -> update positions
					if (lastPhysPos < sourceIt->physicalPosition)
					{//deletion to old reference sequence
						TPos delSize = sourceIt->physicalPosition - lastPhysPos;

						//case 2: deleted entire target node -> skip currentNode
						if (currentNode.length <= delSize)
						{//deleted entire node
							lastPhysPos += currentNode.length;
							shiftVirtualPos -= currentNode.length;
							//access same node in next step - assure the right physical position
							sourceIt->physicalPosition -= shiftPhysicalPos;
							currentNode = *(++targetIt);
							currentNode.virtualPosition += shiftVirtualPos;
							break;
						}
						shiftVirtualPos -= delSize;
						currentNode.physicalPosition += delSize;
						currentNode.length -= delSize;
					}
					//case3: deletion in target to new reference sequence -> right shift of physical position in source
					if (sourceIt->physicalPosition < currentNode.physicalPosition)
					{
						TPos deltaPhysPositions = currentNode.physicalPosition - sourceIt->physicalPosition;
						sourceIt->physicalPosition += deltaPhysPositions;
						shiftPhysicalPos += deltaPhysPositions;
					}
					//case4: cover same area in new reference sequence
					TSize tmpLength = currentNode.length;
					if (tmpLength >  sourceIt->length)		//before node ends a new node in jrn comes before
					{//keep journal node
						currentNode.virtualPosition += sourceIt->length;
						currentNode.physicalPosition += sourceIt->length;
						currentNode.length -= sourceIt->length;
					}
					else
					{
						if (tmpLength < sourceIt->length)  // in reference comes a new node before end of jrn node
						{
							//add new node with new information.
							TJournalEntry tmpNode(*sourceIt);
							tmpNode.virtualPosition = sourceIt->virtualPosition + currentNode.length;
							tmpNode.physicalPosition = sourceIt->physicalPosition + currentNode.length - shiftPhysicalPos; //next node accessed -> assure right phys pos
							tmpNode.length = sourceIt->length - currentNode.length;
							sourceIt->length = currentNode.length;		//change length to new size
							//insert right from old node
							insert((*sourceTreePtr)._journalNodes,pos+1, tmpNode);
							//end of journal nodes changed
							sourceEnd = end(source._journalEntries);
						}
						currentNode = *(++targetIt);
						currentNode.virtualPosition += shiftVirtualPos;
					}
					lastPhysPos = sourceIt->physicalPosition + sourceIt->length;	//last orignal node touched in source
					break;
				}
				case 4u: //tmp->SOURCE_PATCH & jrn->SOURCE_ORIGINAL
				{
					SEQAN_ASSERT_EQ_MSG(currentNode.virtualPosition, sourceIt->virtualPosition,"Source virtual position <%d> does not match target virtual position <%d>!", currentNode.virtualPosition, sourceIt->virtualPosition);
					//update physical position first
					sourceIt->physicalPosition += shiftPhysicalPos;

					//case 1: deleted part of target node -> update positions
					if (lastPhysPos < sourceIt->physicalPosition)
					{//deletion to old reference sequence
						TPos delSize = sourceIt->physicalPosition - lastPhysPos;
						//case 2: deleted entire target node -> skip currentNode
						if (currentNode.length <= delSize)
						{//deleted entire node
							shiftPhysicalPos -= currentNode.length;
							shiftVirtualPos -= currentNode.length;
							sourceIt->physicalPosition -= currentNode.length;
							//access same node in next step - assure the right physical position
							sourceIt->physicalPosition -= shiftPhysicalPos;
							currentNode = *(++targetIt);
							currentNode.virtualPosition += shiftVirtualPos;
							break;
						}
						shiftVirtualPos -= delSize;
						sourceIt->physicalPosition -= delSize;
						shiftPhysicalPos -= delSize;
						currentNode.physicalPosition += delSize;
						currentNode.length -= delSize;
					}

					TSize smallerLength = ( currentNode.length > sourceIt->length) ? sourceIt->length :  currentNode.length;
					shiftPhysicalPos -= smallerLength;
					TInsertionBuffer insBuff = infix(target._insertionBuffer, currentNode.physicalPosition, currentNode.physicalPosition + smallerLength);

					//case4: cover same area in new reference sequence
					TSize tmpLength = currentNode.length;
					if (tmpLength > sourceIt->length)
					{
						currentNode.virtualPosition += smallerLength;
						currentNode.physicalPosition += smallerLength;
						currentNode.length -= sourceIt->length;
					}
					else
					{
						if (tmpLength < sourceIt->length)
						{
							TJournalEntry tmpNode;
							tmpNode.segmentSource = SOURCE_ORIGINAL;
							tmpNode.virtualPosition = sourceIt->virtualPosition + smallerLength;
							tmpNode.physicalPosition = sourceIt->physicalPosition - shiftPhysicalPos;
							tmpNode.length = sourceIt->length - smallerLength;
							insert((*sourceTreePtr)._journalNodes, pos + 1, tmpNode);
							//end of journal nodes changed
							sourceEnd = end(source._journalEntries);
						}
						//take new node of the reference sequence
						currentNode = *(++targetIt);
						currentNode.virtualPosition += shiftVirtualPos;
					}

					sourceIt->segmentSource = SOURCE_PATCH;
					sourceIt->physicalPosition = length(source._insertionBuffer);
					sourceIt->length = smallerLength;
					append(source._insertionBuffer, insBuff, Exact());
					break;
				}
				case 5u: //tmp->SOURCE_ORIGINAL & jrn->SOURCE_PATCH
				{
					SEQAN_ASSERT_EQ_MSG(currentNode.virtualPosition, sourceIt->virtualPosition,"Source virtual position <%d> does not match target virtual position <%d>!", currentNode.virtualPosition, sourceIt->virtualPosition);

					//update shift of virtual position
					shiftVirtualPos += sourceIt->length;
					currentNode.virtualPosition += sourceIt->length;
					break;
				}
				case 6u: //tmp->SOURCE_PATCH & jrn->SOURCE_PATCH
				{
					SEQAN_ASSERT_EQ_MSG(currentNode.virtualPosition, sourceIt->virtualPosition,"Source virtual position <%d> does not match target virtual position <%d>!", currentNode.virtualPosition, sourceIt->virtualPosition);

					//insert before tmp node
					shiftVirtualPos += sourceIt->length;
					currentNode.virtualPosition += sourceIt->length;
					break;
				}
				default:
					SEQAN_ASSERT_FAIL("Failed in %d", __LINE__);
					break;
				}
			}
			else
			{
				++targetIt;
			}
		}

	}

	template <typename TValue, typename THostSpec, typename TBuffSpec>
	void
	synchronizeJournalTrees_(String<TValue, Journaled<THostSpec, UnbalancedTree, TBuffSpec> > const & target,
						  String<TValue, Journaled<THostSpec, UnbalancedTree, TBuffSpec> > & source)
	{
		SEQAN_CHECKPOINT;

		typedef String<TValue, Journaled<THostSpec, UnbalancedTree, TBuffSpec> > TJournalString;
		typedef typename Host<TJournalString>::Type THost;
		typedef String<TValue, TBuffSpec> TInsertionBuffer;

		typedef typename JournalType<TJournalString>::Type TJournalEntries;
		typedef typename TJournalEntries::TNode TNode;
		typedef typename Position<TNode>::Type TPos;
		typedef typename Size<TNode>::Type TSize;
//		typedef JournalEntryLtByVirtualPos<TPos, TSize> TCmp;

		typedef typename Iterator<TJournalEntries, Standard>::Type TIterator;
		typedef typename Iterator<TJournalEntries const, Standard>::Type TConstIterator;

		TConstIterator targetBegin = begin(target._journalEntries);
		TConstIterator targetEnd = end(target._journalEntries);
		TIterator sourceBegin = begin(source._journalEntries);
		TIterator sourceEnd = end(source._journalEntries);
		TIterator targetIt = targetBegin;

		TJournalEntries * sourceTreePtr = &source._journalEntries;
		//manages the current position within the journal entries of the source
		TNode * currentNode = targetBegin._currentNode;

		int shiftPhysicalPos = 0;
		int shiftVirtualPos = 0;
		TPos lastPhysPos = 0;

		//make search different -> search from last accessed node
		TIterator sourceIt;
		TNode * sourceNode;

		//set host information from target to source sequence
		(*sourceTreePtr)._originalStringLength = target._journalEntries._originalStringLength;
		setValue(source._holder, host(target));

		while(targetIt != targetEnd)
		{
			sourceIt = findInJournalEntries((*sourceTreePtr),cargo(*currentNode).virtualPosition);//return const value???
			sourceNode = sourceIt._currentNode;
			//current node points behind last node in the source
			if (cargo(*currentNode).virtualPosition < cargo(*sourceNode).virtualPosition + cargo(*sourceNode).length)
			{
//				TPos pos = sourceIt - sourceBegin;
				switch(cargo(*currentNode).segmentSource + (cargo(*sourceNode).segmentSource << 1))
				{
				case 3u:	//tmp->SOURCE_ORIGINAL & jrn->SOURCE_ORIGINAL
				{
					SEQAN_ASSERT_EQ(cargo(*currentNode).virtualPosition, cargo(*sourceNode).virtualPosition);
					//update physical position first
					cargo(*sourceNode).physicalPosition += shiftPhysicalPos;

					//case 1: deleted part of target node -> update positions
					if (lastPhysPos < cargo(*sourceNode).physicalPosition)
					{//deletion to old reference sequence
						TPos delSize = cargo(*sourceNode).physicalPosition - lastPhysPos;

						//case 2: deleted entire target node -> skip currentNode
						if (cargo(*currentNode).length <= delSize)
						{//deleted entire node
							lastPhysPos += cargo(*currentNode).length;
							shiftVirtualPos -= cargo(*currentNode).length;
							//access same node in next step - assure the right physical position
							cargo(*sourceNode).physicalPosition -= shiftPhysicalPos;
							currentNode = (++targetIt)._currentNode;
							cargo(*currentNode).virtualPosition += shiftVirtualPos;
							break;
						}
						shiftVirtualPos -= delSize;
						cargo(*currentNode).physicalPosition += delSize;
						cargo(*currentNode).length -= delSize;
					}
					//case3: deletion in target to new reference sequence -> right shift of physical position in source
					if (cargo(*sourceNode).physicalPosition < cargo(*currentNode).physicalPosition)
					{
						TPos deltaPhysPositions = cargo(*currentNode).physicalPosition - cargo(*sourceNode).physicalPosition;
						cargo(*sourceNode).physicalPosition += deltaPhysPositions;
						shiftPhysicalPos += deltaPhysPositions;
					}
					//case4: cover same area in new reference sequence
					TSize tmpLength = cargo(*currentNode).length;
					if (tmpLength >  cargo(*sourceNode).length)		//before node ends a new node in jrn comes before
					{//keep journal node
						cargo(*currentNode).virtualPosition += cargo(*sourceNode).length;
						cargo(*currentNode).physicalPosition += cargo(*sourceNode).length;
						cargo(*currentNode).length -= cargo(*sourceNode).length;
					}
					else
					{
						if (tmpLength < cargo(*sourceNode).length)  // in reference comes a new node before end of jrn node
						{
							splitNode(*sourceTreePtr,
									cargo(*sourceNode).virtualPosition + cargo(*currentNode).length,
									cargo(*sourceNode).physicalPosition + cargo(*currentNode).length - shiftPhysicalPos,
									cargo(*sourceNode).length - cargo(*currentNode).length);
							cargo(*sourceNode).length = cargo(*currentNode).length;
							sourceEnd = end(source._journalEntries);
						}
						currentNode = (++targetIt)._currentNode;
						cargo(*currentNode).virtualPosition += shiftVirtualPos;
					}
					lastPhysPos = cargo(*sourceNode).physicalPosition + cargo(*sourceNode).length;	//last orignal node touched in source
					break;
				}
				case 4u: //tmp->SOURCE_PATCH & jrn->SOURCE_ORIGINAL
				{
					SEQAN_ASSERT_EQ(cargo(*currentNode).virtualPosition, cargo(*sourceNode).virtualPosition);
					//update physical position first
					cargo(*sourceNode).physicalPosition += shiftPhysicalPos;

					//case 1: deleted part of target node -> update positions
					if (lastPhysPos < cargo(*sourceNode).physicalPosition)
					{//deletion to old reference sequence
						TPos delSize = cargo(*sourceNode).physicalPosition - lastPhysPos;
						//case 2: deleted entire target node -> skip currentNode
						if (cargo(*currentNode).length <= delSize)
						{//deleted entire node
							shiftPhysicalPos -= cargo(*currentNode).length;
							shiftVirtualPos -= cargo(*currentNode).length;
							cargo(*sourceNode).physicalPosition -= cargo(*currentNode).length;
							//access same node in next step - assure the right physical position
							cargo(*sourceNode).physicalPosition -= shiftPhysicalPos;
							currentNode = (++targetIt)._currentNode;
							cargo(*currentNode).virtualPosition += shiftVirtualPos;
							break;
						}
						shiftVirtualPos -= delSize;
						cargo(*sourceNode).physicalPosition -= delSize;
						shiftPhysicalPos -= delSize;
						cargo(*currentNode).physicalPosition += delSize;
						cargo(*currentNode).length -= delSize;
					}

					TSize smallerLength = ( cargo(*currentNode).length > cargo(*sourceNode).length) ? cargo(*sourceNode).length :  cargo(*currentNode).length;
					shiftPhysicalPos -= smallerLength;
					TInsertionBuffer insBuff = infix(target._insertionBuffer, cargo(*currentNode).physicalPosition, cargo(*currentNode).physicalPosition + smallerLength);

					//case4: cover same area in new reference sequence
					TSize tmpLength = cargo(*currentNode).length;
					if (tmpLength > cargo(*sourceNode).length)
					{
						cargo(*currentNode).virtualPosition += smallerLength;
						cargo(*currentNode).physicalPosition += smallerLength;
						cargo(*currentNode).length -= cargo(*sourceNode).length;
					}
					else
					{
						if (tmpLength < cargo(*sourceNode).length)
						{
							splitNode(*sourceTreePtr,
									cargo(*sourceNode).virtualPosition + smallerLength,
									cargo(*sourceNode).physicalPosition - shiftPhysicalPos,
									cargo(*sourceNode).length - smallerLength);
							sourceEnd = end(source._journalEntries);
						}
						//take new node of the reference sequence
						currentNode = (++targetIt)._currentNode;
						cargo(*currentNode).virtualPosition += shiftVirtualPos;
					}

					cargo(*sourceNode).segmentSource = SOURCE_PATCH;
					cargo(*sourceNode).physicalPosition = length(source._insertionBuffer);
					cargo(*sourceNode).length = smallerLength;
					append(source._insertionBuffer, insBuff, Exact());
					break;
				}
				case 5u: //tmp->SOURCE_ORIGINAL & jrn->SOURCE_PATCH
				{
					SEQAN_ASSERT_EQ(cargo(*currentNode).virtualPosition, cargo(*sourceNode).virtualPosition);
					//update shift of virtual position
					shiftVirtualPos += cargo(*sourceNode).length;
					cargo(*currentNode).virtualPosition += cargo(*sourceNode).length;
					break;
				}
				case 6u: //tmp->SOURCE_PATCH & jrn->SOURCE_PATCH
				{
					SEQAN_ASSERT_EQ(cargo(*currentNode).virtualPosition, cargo(*sourceNode).virtualPosition);
					//insert before tmp node
					shiftVirtualPos += cargo(*sourceNode).length;
					cargo(*currentNode).virtualPosition += cargo(*sourceNode).length;
					break;
				}
				default:
					SEQAN_ASSERT_FAIL("Failed in %d", __LINE__);
					break;
				}
			}
			else
			{
				++targetIt;
			}
		}
	}

	template <typename TCargo>
	inline
	void splitNode(JournalEntries<TCargo, UnbalancedTree> & tree,
				   typename Position<typename JournalEntries<TCargo, UnbalancedTree>::TNode>::Type const & virtualPos,
				   typename Position<typename JournalEntries<TCargo, UnbalancedTree>::TNode>::Type const & physicalBeginPos,
				   typename Size<typename JournalEntries<TCargo, UnbalancedTree>::TNode>::Type const & length)
	{
	    SEQAN_CHECKPOINT;
	    typedef JournalEntries<TCargo, UnbalancedTree> TJournalEntries;
	    typedef typename Iterator<TJournalEntries>::Type TIterator;
	    typedef typename TJournalEntries::TNode TNode;
	    typedef typename Position<TNode>::Type TPos;
	    typedef typename Size<TNode>::Type TSize;

	    SEQAN_ASSERT_TRUE(checkStructure(tree._root));
	    SEQAN_ASSERT_TRUE(checkOrder(tree._root));
	    SEQAN_ASSERT_TRUE(checkVirtualPositions(tree._root));

	    // Handle special case of empty tree.
	    if (tree._root == 0) {
	        SEQAN_ASSERT_EQ(virtualPos, 0u);
	        if (length == 0)
	            return;
	        TNode * tmp;
	        allocate(tree._nodeAllocator, tmp, 1);
	        tree._root = new (tmp) TNode(TCargo(SOURCE_PATCH, physicalBeginPos, virtualPos, length));
	        return;
	    }

	    TNode * node;
	    TNode * parent;
	    TIterator iter = findInJournalEntries(tree, virtualPos);
	    node = iter._currentNode;
	    parent = node->parent;
	    SEQAN_ASSERT_LEQ(cargo(*node).virtualPosition, virtualPos);
	    // Found node that contains virtualPos.
		TPos offset = virtualPos - cargo(*node).virtualPosition;
		TNode * tmp;
		// Current node becomes left part of current node.
		cargo(*node).length = offset;
		// Create insertion node.
		allocate(tree._nodeAllocator, tmp, 1);
		TNode * insertNode = new (tmp) TNode(TCargo(SOURCE_ORIGINAL, physicalBeginPos, virtualPos, length));
		// Insert into tree...
		insertNode->left = node;
		node->parent = insertNode;

		insertNode->right = node->right;
		if (node->right != 0)
			node->right->parent = insertNode;
		node->right = 0;
		if (parent == 0) {
			// current is the root node.
			tree._root = insertNode;
			insertNode->parent = 0;
		} else {
			if (parent->left == node)
				parent->left = insertNode;
			else
				parent->right = insertNode;
			insertNode->parent = parent;
		}

	    SEQAN_ASSERT_TRUE(checkStructure(tree._root));
	    SEQAN_ASSERT_TRUE(checkOrder(tree._root));
	    SEQAN_ASSERT_TRUE(checkVirtualPositions(tree._root));
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TScore, typename TStrategy>
	void
	internalJoin_(String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & source,
				  String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & target,
				  TScore const & /*score*/,
				  TStrategy const & tag)
	{
		SEQAN_CHECKPOINT;

		//TODO rmaerker: extend gap sequence for journal strings
		//TODO rmaerker: allow globalAlignment() function
	}

	template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TScore>
	void
	internalJoin_(String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & source,
				  String <TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & target,
				  TScore const & /*score*/,
				  SimpleJoin const & /*tag*/)
	{
		SEQAN_CHECKPOINT;
		typedef String< TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
		typedef typename JournalType<TJournalString>::Type TJournalType;
		typedef typename Host<TJournalString>::Type THost;
		typename Iterator<TJournalType>::Type it = begin(source._journalEntries);

		//temporary sequence of the journal
		THost tmpHost;
		for (;it != end(source._journalEntries);++it)
		{
			if (value(it).segmentSource == SOURCE_ORIGINAL)
			{
				append(tmpHost,infix(value(source._holder), value(it).physicalPosition, value(it).physicalPosition + value(it).length));
				continue;
			}
			SEQAN_ASSERT_EQ(value(it).segmentSource,SOURCE_PATCH);
			append(tmpHost,infix(source._insertionBuffer, value(it).physicalPosition, value(it).physicalPosition + value(it).length));
		}

		//set the host of the group reference sequence to the new host of the given journal
		// and insert the temporary string at position 0 to the journal.
		// Delete the host sequence.
		setHost(source, host(target));
		insert(source,0,tmpHost);
		erase(source,length(tmpHost),length(source));
	}



	template <typename TScoreValue, unsigned DIMENSION, typename TString>
	TScoreValue
	computeSmallestJournal_(Matrix<TScoreValue, DIMENSION> & diag_matrix_,
		Matrix<TScoreValue, DIMENSION> & vert_matrix_,
		Matrix<TScoreValue, DIMENSION> & hori_matrix_,
		TString const & target_,
		TString const & source_,
		Score<TScoreValue, ExtendedScore> const & score_)
	{
		SEQAN_CHECKPOINT;

		typedef Matrix<TScoreValue, DIMENSION> TMatrix;

		typedef typename Size<TMatrix>::Type TSize;
		typedef typename Iterator<TMatrix, Rooted>::Type TMatrixIterator;


		typedef typename Iterator<TString, Standard>::Type TStringIterator;
		typedef typename Position<TString>::Type TPosition;
		typedef typename Value<TString const>::Type TValue;

		//-------------------------------------------------------------------------
		//define some variables
		TSize target_length = length(target_);
		TSize source_length = length(source_);

		TPosition x_begin = 0;
		TPosition x_end = length(target_);
		TPosition y_begin = 0;
		TPosition y_end = length(source_);
		TPosition const inf = 100000;

		//define scores for recursion
		TScoreValue score_match = scoreMatch(score_);
		TScoreValue score_mismatch = scoreMisMatch(score_);
		TScoreValue score_gap_open_insert = scoreGapOpenInsert(score_);
		TScoreValue score_gap_extend_insert = scoreGapExtendInsert(score_);
		TScoreValue score_gap_open_delete = scoreGapOpenDelete(score_);
		TScoreValue score_gap_extend_delete = scoreGapExtendDelete(score_);

		TScoreValue borderHori_ = score_gap_open_insert;
		TScoreValue borderHoriExtend_ = score_gap_extend_insert;
		TScoreValue borderVert_ = score_gap_open_delete;
		TScoreValue borderVertExtend_ = score_gap_extend_delete;
		TScoreValue v;

		setDimension(diag_matrix_, 2);
		setLength(diag_matrix_, 0, target_length + 1);	//x-coordinate
		setLength(diag_matrix_, 1, source_length + 1);	//y-coordinate
		resize(diag_matrix_);
		setDimension(vert_matrix_, 2);
		setLength(vert_matrix_, 0, target_length + 1);	//x-coordinate
		setLength(vert_matrix_, 1, source_length + 1);	//y-coordinate
		resize(vert_matrix_);
		setDimension(hori_matrix_, 2);
		setLength(hori_matrix_, 0, target_length + 1);	//x-coordinate
		setLength(hori_matrix_, 1, source_length + 1);	//y-coordinate
		resize(hori_matrix_);

		TMatrixIterator diag_col_ = end(diag_matrix_) - 1;
		TMatrixIterator diag_finger_1_;
		TMatrixIterator diag_finger_2_;
		TMatrixIterator vert_col_ = end(vert_matrix_) - 1;
		TMatrixIterator vert_finger_1_;
		TMatrixIterator vert_finger_2_;
		TMatrixIterator hori_col_ = end(hori_matrix_) - 1;
		TMatrixIterator hori_finger_1_;
		TMatrixIterator hori_finger_2_;

		//-------------------------------------------------------------------------
		// init

		diag_finger_1_ = diag_col_;
		*diag_finger_1_ = 0;
		vert_finger_1_ = vert_col_;
		*vert_finger_1_ = inf;
		hori_finger_1_ = hori_col_;
		*hori_finger_1_ = inf;
		for (TPosition x = x_end; x > x_begin; --x)
		{
			goPrevious(diag_finger_1_, 0);
			*diag_finger_1_ = borderVert_;
			goPrevious(hori_finger_1_, 0);
			*hori_finger_1_ = inf;
			goPrevious(vert_finger_1_, 0);
			*vert_finger_1_ = borderVert_;
			borderVert_ += borderVertExtend_;
		}
		//------------------------------------------------------------------------
		//fill matrix
		//borderVert_ = score_gap_open_insert;
		for (TPosition y = y_end; y > y_begin; --y)
		{
			TValue cy = source_[y-1];
			v = borderHori_;

			vert_finger_2_ = vert_col_;	//points to last column
			goPrevious(vert_col_, 1);	//points to this column
			vert_finger_1_ = vert_col_;
			*vert_finger_1_ = inf;          //initialize first column

			diag_finger_2_ = diag_col_;
			goPrevious(diag_col_, 1);
			diag_finger_1_ = diag_col_;
			*diag_finger_1_ = v;			//initialize first column

			hori_finger_2_ = hori_col_;
			goPrevious(hori_col_, 1);
			hori_finger_1_ = hori_col_;
			*hori_finger_1_ = v;	//initialize first column

			for (TPosition x = x_end; x > x_begin; --x)
			{//allow no mismatches
				//compute entry in diag_matrix
				goPrevious(diag_finger_1_, 0);
				v = (*diag_finger_2_ < *vert_finger_2_) ? *diag_finger_2_ : *vert_finger_2_;
				v = (*hori_finger_2_ < v) ? *hori_finger_2_ : v;
				if (target_[x - 1] == cy) *diag_finger_1_ = v + score_match;
				else *diag_finger_1_ = v + score_mismatch;

				//compute entry in vert_matrix
				v = *vert_finger_1_;
				goPrevious(vert_finger_1_, 0);
				v += score_gap_extend_delete;
				goPrevious(diag_finger_2_, 1);
				*vert_finger_1_ = (v < (*diag_finger_2_ + score_gap_open_delete)) ? v : (*diag_finger_2_ + score_gap_open_delete);
				//what happens with the horizontal matrix?
				goPrevious(vert_finger_2_, 0);

				//compute entry in hori_matrix
				goPrevious(hori_finger_2_, 0);
				goPrevious(hori_finger_1_, 0);
				v = *hori_finger_2_ + score_gap_extend_insert;
				goNext(diag_finger_2_, 1);
				goPrevious(diag_finger_2_, 0);
				*hori_finger_1_ = (v < (*diag_finger_2_ + score_gap_open_insert)) ? v : (*diag_finger_2_ + score_gap_open_insert);
				//what happens with the horizontal matrix?
			}
			borderHori_ += borderHoriExtend_;
		}

		v = (*vert_finger_1_ < *hori_finger_1_) ? *vert_finger_1_ : *hori_finger_1_;
		v = (*diag_finger_1_ < v) ? *diag_finger_1_ : v;

		return v;
	}

	template <typename TScoreValue, unsigned DIMENSION,  typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec, typename TInfoNodes>
	void
	traceSmallestJournal_(Matrix<TScoreValue, DIMENSION>  & diagMatrix,
				  Matrix<TScoreValue, DIMENSION>  & vertMatrix,
				  Matrix<TScoreValue, DIMENSION>  & horiMatrix,
				  String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & target,
				  String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > const & source,
				  TInfoNodes & tmpNodes,
				  Score<TScoreValue, ExtendedScore> const & score)
	{
		SEQAN_CHECKPOINT;

		typedef Iter<Matrix<TScoreValue, DIMENSION>, PositionIterator > TMatrixIterator;
		typedef typename Position<Matrix<TScoreValue, DIMENSION> >::Type TPosition;

		typedef String<TValue, Journaled<THostSpec,TJournalSpec,TBuffSpec> > TJournalString;
		typedef typename Position<TJournalString>::Type TJournalPosition;

		typedef typename Iterator<TJournalString const, Standard>::Type TJournalIterator;
		typedef typename Host<TJournalString>::Type THost;

		TJournalIterator it_target = begin(target);
		TJournalIterator it_target_end = end(target);
		TJournalIterator it_source = begin(source);
		TJournalIterator it_source_end = end(source);

		TScoreValue score_gap_open_insert = scoreGapOpenInsert(score);
		TScoreValue score_gap_open_delete = scoreGapOpenDelete(score);

		TMatrixIterator diag_ = begin(diagMatrix);	//points to the first entry in the diagonal matrix
		TMatrixIterator hori_ = begin(horiMatrix);	//points to the first entry in the horizontal matrix
		TMatrixIterator vert_ = begin(vertMatrix);	//points to the first entry in the vertival matrix
		TPosition pos = position(diag_);

		// indicate which matrix we are in
		bool hori = false, vert = false, diag = false;
		if (*diag_ < *hori_)
		{
			if (*diag_ < *vert_) diag = true;
			else vert = true;
		}
		else
		{
			if (*hori_ < *vert_) hori = true;
			else vert = true;
		}

		TJournalPosition delLength = 0;	//stores the length of the deletion
		int delta = 0;
		THost insertionBuff;		//stores the insertion chars note it is reversed since the trace back is started at the end
		//-------------------------------------------------------------------------
		//follow the trace until the border is reached
		while (pos != position(end(diagMatrix))-1)
		{
			if(diag)
			{//is in diagonal matrix
				++it_target;
				++it_source;
//				pos += dim_0_len + 1; //position in the matrix
				goNext(diag_,0);
				goNext(diag_,1);
				pos = position(diag_);
				if (getValue(diagMatrix,pos) <= getValue(horiMatrix,pos))
				{//d < h
					if (getValue(vertMatrix,pos) < getValue(diagMatrix,pos))
					{ // v < d < h
						vert = true;
						diag = false;
					}
				}
				else
				{//h < d
					diag = false;
					hori = true;
					if (getValue(vertMatrix,pos) < getValue(horiMatrix,pos))
					{//h < v < d
						hori = false;
						vert = true;
					}
				}
			}
			else
			{//is not in diagonal matrix
				if(hori)
				{// is in horizontal matrix -> recordInsertion
					appendValue(insertionBuff, value(it_source));
					++it_source;
//					pos += dim_0_len;
					goNext(diag_,1);
					pos = position(diag_);

					if ( getValue(vertMatrix,pos) + score_gap_open_insert < getValue(horiMatrix,pos)) // continue in vertical matrix
					{
						vert = true;
						hori = false;

						TJournalPosition insPos = coordinate(diag_,0);
//						insert(source2,insPos + delta, insertionBuff);
						recordInsertion_(tmpNodes, insPos, delta, insertionBuff);

					}
					else if ( getValue(diagMatrix,pos) + score_gap_open_insert < getValue(horiMatrix,pos)) // continue in diagonal matrix
					{
						diag = true;
						hori = false;

						TJournalPosition insPos = coordinate(diag_,0);
//						insert(source2,insPos + delta, insertionBuff);
						recordInsertion_(tmpNodes, insPos, delta, insertionBuff);
					}
				}
				else
				{//is in vertical matrix -> record deletion
					if(vert)
					{//check is in horizontal
						++delLength;
						++it_target;
//						++pos;
						goNext(diag_,0);
						pos = position(diag_);
						if (getValue(horiMatrix,pos) + score_gap_open_delete < getValue(vertMatrix,pos))
						{
							hori = true;
							vert = false;

							TJournalPosition delPos = coordinate(diag_,0);
							recordDeletion_(tmpNodes,delPos, delLength, delta, insertionBuff);
//							erase(source2, delPos - delLength, delPos);
//							delta -= delLength;
//							delLength = 0;
						}
						else if (getValue(diagMatrix,pos) + score_gap_open_delete < getValue(vertMatrix,pos) )
						{
							diag = true;
							vert = false;

							TJournalPosition delPos = coordinate(diag_,0);
							recordDeletion_(tmpNodes,delPos, delLength, delta,insertionBuff);
//							erase(source2, delPos - delLength, delPos);
//							delta -= delLength;
//							delLength = 0;
						}
					}
				}
			}
		}
	}

	template <typename TInfoNodes, typename TPos, typename TPos2, typename TInsBuff>
	inline void
	recordInsertion_(TInfoNodes & nodes,
					 TPos & beginPos,
					 TPos2 & delta,
					 TInsBuff & insBuff)
	{
		SEQAN_CHECKPOINT;

		typedef typename Value<TInfoNodes>::Type TInfoNode;

		TInfoNode info;
		info.nodeType_ = NODE_INSERTION;
		info.begin_ = beginPos + delta;
		info.length_ = length(insBuff);
		info.insertionBuff_ = insBuff;
		appendValue(nodes, info, Exact());

		delta = length(insBuff);
		clear(insBuff);
	}

	template <typename TInfoNodes, typename TPos, typename TPos2, typename TInsBuff>
	inline void
	recordDeletion_(TInfoNodes & nodes,
					 TPos & endPos,
					 TPos & delLength,
					 TPos2 & delta,
					 TInsBuff & /*insBuff*/)
	{
		SEQAN_CHECKPOINT;

		typedef typename Value<TInfoNodes>::Type TInfoNode;

		TInfoNode info;
		info.nodeType_ = NODE_DELETION;
		info.begin_ = endPos - delLength + delta;
		info.length_ = delLength;
		appendValue(nodes, info, Exact());

		delta -= delLength;
		delLength = 0;
	}
}

#endif /* STRING_SET_JOURNALED_UTIL_H_ */
