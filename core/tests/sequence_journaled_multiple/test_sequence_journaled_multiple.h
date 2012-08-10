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
  Defines test for the sequence_journaled_multiple package
  ==========================================================================*/

#ifndef TEST_SEQUENCE_JOURNALED_MULTIPLE_H_
#define TEST_SEQUENCE_JOURNALED_MULTIPLE_H_

#include <cstdlib>
#include <sstream>
#include <string>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journaled.h>
#include <seqan/sequence_journaled_multiple.h>

using namespace seqan;

template <typename TSet, typename TSpec>
void testStringSetJournaledGroup()
{
	typedef typename Position<TSet>::Type TPosition;
	Group<TSet, TSpec> group1;
	typedef typename GroupLabel<Group<TSet, TSpec> >::Type TLabel;
	typedef typename GroupMembers<Group<TSet, TSpec> >::Type TMembers;

	TPosition pos1 = 0;
	TPosition pos2 = 1;
	TPosition pos3 = 2;
	TPosition pos4 = 3;
	TPosition pos5 = 4;
	TPosition pos6 = 5;

	TLabel lbl1 = "test1";
	TMembers mem1;
	resize(mem1, 4);
	mem1[0] = pos1;
	mem1[1] = pos2;
	mem1[2] = pos3;
	mem1[3] = pos4;

	setGroupLabel(group1,lbl1);
	setGroupMembers(group1, mem1);
	setGroupReference(group1, pos2);

	SEQAN_ASSERT_EQ(getGroupLabel(group1),lbl1);
	SEQAN_ASSERT_EQ(length(group1),4);
	SEQAN_ASSERT_EQ(getGroupReference(group1),pos2);
	SEQAN_ASSERT_TRUE(group1[0] == pos1);


	SEQAN_ASSERT_TRUE(join(group1, pos5));
	SEQAN_ASSERT_TRUE(isMember(group1, pos5));
	SEQAN_ASSERT_NOT(join(group1, pos3));

	SEQAN_ASSERT_TRUE(join(group1, pos6));
	SEQAN_ASSERT_TRUE(isMember(group1, pos6));
	SEQAN_ASSERT_EQ(length(group1), 6);

	SEQAN_ASSERT_NOT(swap(group1, 10u));
	SEQAN_ASSERT_TRUE(swap(group1, pos4));
	SEQAN_ASSERT_TRUE(isReference(group1,pos4));
	SEQAN_ASSERT_EQ(length(group1),6);

	SEQAN_ASSERT_NOT(disjoin(group1, (TPosition) 10));
	SEQAN_ASSERT_TRUE(disjoin(group1, pos1));
	SEQAN_ASSERT_EQ(length(group1), 5);
	SEQAN_ASSERT_NOT(disjoin(group1, getGroupReference(group1)));

	Group<TSet, TSpec> group2;
	group2 = group1;

	SEQAN_ASSERT_TRUE(group2 == group1);

	TLabel lbl2 = "test2";

	rename(group2,lbl2);

	SEQAN_ASSERT_EQ(getGroupLabel(group2), lbl2);
	SEQAN_ASSERT_TRUE(group2 != group1);

	Group<TSet, TSpec> const group3(group1);
	SEQAN_ASSERT_TRUE(group3 != group2);
	SEQAN_ASSERT_TRUE(group3 == group1);
	SEQAN_ASSERT_TRUE(group3[0] == pos2);
	SEQAN_ASSERT_EQ(getGroupLabel(group3),lbl1);
	SEQAN_ASSERT_EQ(length(group3),5);
	SEQAN_ASSERT_EQ(getGroupReference(group3),pos4);

	clear(group2);
	SEQAN_ASSERT_EQ(getGroupReference(group2), MaxValue<TPosition>::VALUE);
	SEQAN_ASSERT_TRUE(length(group2) == 0);
	SEQAN_ASSERT_TRUE(getGroupLabel(group2) == "");
}

template <typename TValue, typename TSetSpec>
void
testStringSetBasic()
{
	//TODO rmaerker: define tests for the  string set
	//test appendValue 1 und 2
	//test resize
	//test resolveGroup
	//test length
	//test copy-ctor
	//test standard ctor
	//test assignment operator
	//test getter and setter functions for members
	//test assignValue function

	typedef typename Host<TValue>::Type THost;

	typedef StringSet<TValue, Owner<TSetSpec> > TJournalSet;

	TJournalSet jrnSet1;

	SEQAN_ASSERT_EQ(length(jrnSet1), 0);
	SEQAN_ASSERT_EQ(length(getStrings(jrnSet1)), 0);
	SEQAN_ASSERT_EQ(length(getGroups(jrnSet1)), 0);
	SEQAN_ASSERT_EQ(length(getGroupMap(jrnSet1)), 0);

	THost host1 = "host1";
	TValue jrn1(host1);

	appendValue(jrnSet1, jrn1);

	SEQAN_ASSERT_EQ(length(jrnSet1), 1u);
	SEQAN_ASSERT_EQ(length(getStrings(jrnSet1)), 1u);
	SEQAN_ASSERT_EQ(length(getGroups(jrnSet1)), 1u);
	SEQAN_ASSERT_EQ(getGroupReference(getGroups(jrnSet1)[0]),0u);
	SEQAN_ASSERT_EQ(length(getGroupMap(jrnSet1)), 1u);
	SEQAN_ASSERT_EQ(getGroupMap(jrnSet1)[0],0u);

	TValue jrn2(host1);
	insert(jrn2, 2, "insert");

	appendValue(jrnSet1, jrn2);

	SEQAN_ASSERT_EQ(length(jrnSet1), 2u);
	SEQAN_ASSERT_EQ(length(getStrings(jrnSet1)), 2u);
	SEQAN_ASSERT_EQ(length(getGroups(jrnSet1)), 2u);
	SEQAN_ASSERT_EQ(getGroupReference(getGroups(jrnSet1)[1]),1u);
	SEQAN_ASSERT_EQ(length(getGroupMap(jrnSet1)), 2u);
	SEQAN_ASSERT_EQ(getGroupMap(jrnSet1)[1],1u);

	THost host2 = "rost";
	appendValue(jrnSet1, host2);
	SEQAN_ASSERT_EQ(length(jrnSet1), 3u);
	SEQAN_ASSERT_EQ(length(getStrings(jrnSet1)), 3u);
	SEQAN_ASSERT_EQ(length(getGroups(jrnSet1)), 3u);
	SEQAN_ASSERT_EQ(getGroupReference(getGroups(jrnSet1)[2]),2u);
	SEQAN_ASSERT_EQ(length(getGroupMap(jrnSet1)), 3u);
	SEQAN_ASSERT_EQ(getGroupMap(jrnSet1)[2],2u);


}

template <typename TSource>
void
testSequenceJournaledIsEmpty()
{
	typedef typename Host<TSource>::Type THost;
	typedef typename Position<TSource>::Type TPosition;

	THost src = "empty journal";
	TSource emptyJrn(src);

	SEQAN_ASSERT_TRUE(isFlat(emptyJrn));

	THost ins = "is not an";
	insert(emptyJrn,0, ins);

	SEQAN_ASSERT_NOT(isFlat(emptyJrn));

	TPosition endP = length(ins);
	erase(emptyJrn, (TPosition)0, endP);

	SEQAN_ASSERT_TRUE(isFlat(emptyJrn));

	erase(emptyJrn,0);
	SEQAN_ASSERT_NOT(isFlat(emptyJrn));
}

template <typename TJournal>
void test()
{
	typedef typename Host<TJournal>::Type THost;

	THost host1 = "ABCD";
	TJournal jrn1(host1);
	TJournal jrn2 = jrn1;
	THost ins1 = "XXX";
	insert(jrn1, 2, ins1);
	std::cout <<"JournalEntries: " << jrn1._journalEntries << std::endl;
	THost ins2 = "YY";
	insert(jrn1, 3, ins2);
	std::cout <<"JournalEntries: " << jrn1._journalEntries << std::endl;

	erase(jrn1, 2, 7);
	std::cout <<"JournalEntries: " << jrn1._journalEntries << std::endl;

	std::cout <<"JournalEntries: " << jrn2._journalEntries << std::endl;
	erase(jrn2, 1, 3);
	std::cout <<"JournalEntries: " << jrn2._journalEntries << std::endl;
}

template <typename TJournalString>
inline void
testStringSetJournaledJoinMaxCompression()
{
	typedef typename Host<TJournalString>::Type THost;
	typedef typename Position<TJournalString>::Type TPosition;

	THost sourceHost = "AAACCAAAAA";
	THost targetHost = "AAAAATTAAAAA";

	TJournalString source(sourceHost);
	TJournalString target(targetHost);

	Score<long,ExtendedScore> score(0,-100000,-25,-1,-12,0);
	internalJoin_(source, target, score, SimpleJoin());
	std::cout << source << " " << host(source)<<std::endl;
}

template <typename TJournalString>
void
testStringSetJournaledMerge()
{
	typedef typename Host<TJournalString>::Type THost;
	typedef typename Position<TJournalString>::Type TPosition;

	//case 1: deletion at beginning, the end and in the middle of the old reference sequence
	//journal is same as old reference
	THost newRef1 = "BBBAABBBBAAABBBBB";
	THost oldRef1 = "AAAAA";

	TJournalString newReference1(newRef1);
	TJournalString member1(oldRef1);
	erase(newReference1, 0, 3);
	erase(newReference1, 2, 6);
	erase(newReference1, 5, 10);

	synchronizeJournalTrees_(newReference1, member1);

	SEQAN_ASSERT_TRUE(member1 == "AAAAA");
	SEQAN_ASSERT_TRUE(host(member1) == newRef1);

	//case 2: deletion at beginning, the end and in the middle of the journal sequence
	//old reference is same as new reference

	THost newRef2 = "BBBAABBBBAAABBBBB";
	THost oldRef2 = "BBBAABBBBAAABBBBB";

	TJournalString newReference2(newRef2);
	TJournalString member2(oldRef2);
	erase(member2, 0, 3);
	erase(member2, 2, 6);
	erase(member2, 5, 10);

	synchronizeJournalTrees_(newReference2, member2);

	SEQAN_ASSERT_TRUE(member2 == "AAAAA");
	SEQAN_ASSERT_TRUE(host(member2) == newRef2);

	//case 3: insertion at beginning, the end and in the middle of the reference sequence
	//journal is same as old reference

	THost newRef3 = "AAAAA";
	THost oldRef3 = "BBBAABBBBAAABBBBB";

	TJournalString newReference3(newRef3);
	TJournalString member3(oldRef3);
	insert(newReference3, 0, "BBB");
	insert(newReference3, 5, "BBBB");
	insert(newReference3, 12, "BBBBB");

	synchronizeJournalTrees_(newReference3, member3);

	SEQAN_ASSERT_TRUE(member3 == "BBBAABBBBAAABBBBB");
	SEQAN_ASSERT_TRUE(host(member3) == newRef3);

	//case 4: insertion at beginning, the end and in the middle of the journal sequence
	//new reference is same as old reference

	THost newRef4 = "AAAAA";
	THost oldRef4 = "AAAAA";

	TJournalString newReference4(newRef4);
	TJournalString member4(oldRef4);
	insert(member4, 0, "BBB");
	insert(member4, 5, "BBBB");
	insert(member4, 12, "BBBBB");

	synchronizeJournalTrees_(newReference4, member4);

	SEQAN_ASSERT_TRUE(member4 == "BBBAABBBBAAABBBBB");
	SEQAN_ASSERT_TRUE(host(member4) == newRef4);

	//case 5: deletion from old reference to new reference
	//deletion from journal to old reference

	THost newRefHost5 = "XXXBBBAAXXBBBBXXAAAXBB";
	THost oldRefHost5 = "XXXAAXXXXAAAX";

	TJournalString oldReference5(newRefHost5); //XXXAAXXXXAAAX
	erase(oldReference5,20,22);
	erase(oldReference5,10,14);
	erase(oldReference5,3,6);

	TJournalString member5(oldRefHost5);	//XAAXXAX
	erase(member5,11,13);
	erase(member5,5,7);
	erase(member5,0,2);

	synchronizeJournalTrees_(oldReference5, member5);

	SEQAN_ASSERT_TRUE(member5 == "XAAXXAA");
	SEQAN_ASSERT_TRUE(host(member5) == newRefHost5);

	//case 6: insertion from old reference to new reference
	//insertion from journal to old reference

	THost newRefHost6 = "AAAAA";
	THost oldRefHost6 = "XXXAAXXXXAAAX";

	TJournalString oldReference6(newRefHost6); //XXXAAXXXXAAAX
	insert(oldReference6,0,"XXX");
	insert(oldReference6,5,"XXXX");
	insert(oldReference6,12,"X");

	TJournalString member6(oldRefHost6);	//"XXXBBBAAXXBBBBXXAAAXBB"
	insert(member6,3, "BBB");
	insert(member6,10,"BBBB");
	insert(member6,20,"BB");

	synchronizeJournalTrees_(oldReference6, member6);

	SEQAN_ASSERT_TRUE(member6 == "XXXBBBAAXXBBBBXXAAAXBB");
	SEQAN_ASSERT_TRUE(host(member6) == newRefHost6);

	THost newRefHost7 = "OODDDDOOOOODDDDOODDD";
	THost oldRefHost7 = "1234OOOOOOO1256734OO";

	TJournalString oldReference7(newRefHost7); //1234OOOOOOO1256734OO
	insert(oldReference7,0,"1234");
	erase(oldReference7,6,10);
	insert(oldReference7,11,"1234");
	insert(oldReference7, 13, "567");
	erase(oldReference7, 18, 22);
	erase(oldReference7,20,23);
											//"--123-4OO-----OOOOO12567------34OO";
											//"0123456789012345678901234567890123
	TJournalString member7(oldRefHost7);	//"JJ--3J---JJJJJ-------567JJJJJJ34O-"
	insert(member7,0, "JJ");
	erase(member7,2,4);
	insert(member7,3, "J");
	erase(member7,4, 7);
	insert(member7,4,"JJJJJ");
	erase(member7,9, 16);
	insert(member7,12,"JJJJJJ");
	erase(member7,21,22);

	SEQAN_ASSERT_TRUE(member7 == "JJ3JJJJJJ567JJJJJJ34O");

	synchronizeJournalTrees_(oldReference7, member7);

	SEQAN_ASSERT_TRUE(member7 == "JJ3JJJJJJ567JJJJJJ34O");
	SEQAN_ASSERT_TRUE(host(member7) == newRefHost7);

	THost newRefHost8 = "OODDDDOOOOODDDDOODDD";
	THost oldRefHost8 = "IIIIOOOOOOOIIIIIIIOO";

	TJournalString oldReference8(newRefHost8); //IIIIOOOOOOOIIIIIIIOO
	insert(oldReference8,0,"IIII");
	erase(oldReference8,6,10);
	insert(oldReference8,11,"IIII");
	insert(oldReference8, 13, "III");
	erase(oldReference8, 18, 22);
	erase(oldReference8,20,23);

	TJournalString member8(oldRefHost8);	//"JJIJJJJJJIIIJJJJJJIIO"
	insert(member8,0, "JJIJJJJJJIIIJJJJJJIIO");
	erase(member8,21,41);

	SEQAN_ASSERT_TRUE(member8 == "JJIJJJJJJIIIJJJJJJIIO");

	synchronizeJournalTrees_(oldReference8, member8);

	SEQAN_ASSERT_TRUE(member8 == "JJIJJJJJJIIIJJJJJJIIO");
	SEQAN_ASSERT_TRUE(host(member8) == newRefHost8);

	THost newRefHost9 = "OODDDDOOOOODDDDOODDD";
	THost oldRefHost9 = "1234OOOOOOO1256734OO";

	TJournalString oldReference9(newRefHost9); //1234OOOOOOO1256734OO
	insert(oldReference9,0,"1234OOOOOOO1256734OO");
	erase(oldReference9,20,40);
											//"--123-4OO-----OOOOO12567------34OO";
											//"0123456789012345678901234567890123
	TJournalString member9(oldRefHost9);	//"JJ--3J---JJJJJ-------567JJJJJJ34O-"
	insert(member9,0, "JJ");
	erase(member9,2,4);
	insert(member9,3, "J");
	erase(member9,4, 7);
	insert(member9,4,"JJJJJ");
	erase(member9,9, 16);
	insert(member9,12,"JJJJJJ");
	erase(member9,21,22);

	SEQAN_ASSERT_TRUE(member9 == "JJ3JJJJJJ567JJJJJJ34O");

	synchronizeJournalTrees_(oldReference9, member9);

	SEQAN_ASSERT_TRUE(member9 == "JJ3JJJJJJ567JJJJJJ34O");
	SEQAN_ASSERT_TRUE(host(member9) == newRefHost9);
}

SEQAN_DEFINE_TEST(test_string_set_journaled_group_basic)
{
	testStringSetJournaledGroup< String<size_t>, GroupBasic >();
}

SEQAN_DEFINE_TEST(test_string_set_journaled_basic)
{
	testStringSetBasic< String< char, Journaled<Alloc<void>, SortedArray, Alloc<void> > >, JournalSetBasic>();
}

SEQAN_DEFINE_TEST(test_sequence_journaled_is_empty)
{
	testSequenceJournaledIsEmpty< String<char, Journaled< Alloc<void>, SortedArray > > >();
}

SEQAN_DEFINE_TEST(test_string_set_journaled_uitl)
{
//	testStringSetJournaledUtil<String<char, Journaled<Alloc<void>, SortedArray > > >();
	testStringSetJournaledMerge<String<char, Journaled<Alloc<void>, SortedArray> > >();
	testStringSetJournaledMerge<String<char, Journaled<Alloc<void>, UnbalancedTree> > >();
}

SEQAN_DEFINE_TEST(test_string_set_journaled_join)
{
	testStringSetJournaledJoinMaxCompression<String<char, Journaled<Alloc<void>, SortedArray > > >();
}

#endif /* TEST_SEQUENCE_JOURNALED_MULTIPLE_H_ */
