// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef TESTS_ALIGN_TEST_ALIGN_GAPS_H_
#define TESTS_ALIGN_TEST_ALIGN_GAPS_H_

#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_TEST

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/basic.h>


using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSpec>
void TestGapsBase()
{
	typedef Gaps<TSource, TSpec> TGaps;

	TGaps gaps1;					//default ctor
	TSource src1 = "hello";
	setSource(gaps1, src1);			//setSource
    SEQAN_ASSERT(source(gaps1) == src1);
    SEQAN_ASSERT(getObjectId(source(gaps1)) == getObjectId(src1));

    SEQAN_ASSERT(getObjectId(source(gaps1)) == getObjectId(gaps1));   //id;

	assignSource(gaps1, "blabla");	//assignSource
    SEQAN_ASSERT(source(gaps1) == "blabla");
    SEQAN_ASSERT(src1 == "blabla");

	assignSource(gaps1, "abcdef");	//assignSource
	setClippedBeginPosition(gaps1, 1);
	setBeginPosition(gaps1, 0);
	setClippedEndPosition(gaps1, 5);
    SEQAN_ASSERT(source(gaps1) == "abcdef");
	SEQAN_ASSERT(sourceSegment(gaps1) == "bcde"); //sourceSegment

	moveSource(gaps1, "hullahulla");	//moveSource
    SEQAN_ASSERT(source(gaps1) == "hullahulla");

	moveSource(gaps1, "abcdef", 1, 5);	//moveSource
	SEQAN_ASSERT(source(gaps1) == "abcdef"); //???Sollte das anders sein?
    SEQAN_ASSERT(sourceSegment(gaps1) == "bcde");

	detach(gaps1);			//detach, createSource
    SEQAN_ASSERT(source(gaps1) == src1);
    SEQAN_ASSERT(getObjectId(gaps1) != getObjectId(src1));

	TGaps gaps2(gaps1);				//copy ctor

//  it is a real copy
//	SEQAN_ASSERT(getObjectId(gaps1) == getObjectId(gaps2)) //(its not a real copy)

//____________________________________________________________________________

	setSource(gaps1, src1);
	src1 = "hello";
    SEQAN_ASSERT(getObjectId(gaps1) != getObjectId(gaps2));
	gaps2 = gaps1;					//operator =
    SEQAN_ASSERT(getObjectId(gaps1) == getObjectId(gaps2));
    SEQAN_ASSERT(getObjectId(gaps2) == getObjectId(src1));

	TGaps gaps3(src1);				//ctor with source
	TGaps const & c_gaps3 = gaps3;	//const version

    SEQAN_ASSERT(getObjectId(gaps3) == getObjectId(src1));
    SEQAN_ASSERT(getObjectId(c_gaps3) == getObjectId(src1));

	SEQAN_ASSERT(dependentSource(gaps3));	//dependentSource
    SEQAN_ASSERT(dependentSource(c_gaps3));

	SEQAN_ASSERT(length(gaps3) == length(src1));		//length
	SEQAN_ASSERT(sourceLength(gaps3) == length(src1));	//sourceLength

	SEQAN_ASSERT(clippedBeginPosition(gaps3) == 0);		//clippedBeginPosition
	SEQAN_ASSERT(sourceBegin(gaps3) == begin(source(gaps3)));		//sourceBegin
	SEQAN_ASSERT(sourceBegin(gaps3, Rooted()) == begin(source(gaps3), Rooted()));

	SEQAN_ASSERT(clippedEndPosition(gaps3) == length(src1));	//clippedEndPosition
	SEQAN_ASSERT(sourceEnd(gaps3) == end(source(gaps3)));		//sourceEnd
	SEQAN_ASSERT(sourceEnd(gaps3, Rooted()) == end(source(gaps3), Rooted()));

	SEQAN_ASSERT(*(--end(gaps3)) == 'o'); //end
	SEQAN_ASSERT(*(--end(c_gaps3)) == 'o'); //end

	setClippedBeginPosition(gaps3, 3); //"---lo"			//setClippedBeginPosition
	SEQAN_ASSERT(gaps3 == "lo");
	SEQAN_ASSERT_EQ(clippedBeginPosition(gaps3), 3u);		//clippedBeginPosition
	SEQAN_ASSERT_EQ(beginPosition(gaps3), 3u);			//beginPosition
	SEQAN_ASSERT_EQ(length(gaps3), 2u);					//length
	
	setClippedBeginPosition(gaps3, 1); //"-ello"			//setClippedBeginPosition
	SEQAN_ASSERT(gaps3 == "ello");
	SEQAN_ASSERT_EQ(clippedBeginPosition(gaps3), 1u);		//clippedBeginPosition
	SEQAN_ASSERT_EQ(beginPosition(gaps3), 1u);			//beginPosition
	SEQAN_ASSERT_EQ(length(gaps3),  4u);					//length
	SEQAN_ASSERT(isGap(gaps3, 0u));						//isGap
    SEQAN_ASSERT(gaps3[1] == 'e');
    SEQAN_ASSERT(c_gaps3[1] == 'e');

	SEQAN_ASSERT(getValue(gaps3, 1) == 'e');			//getValue
	SEQAN_ASSERT(getValue(c_gaps3, 1) == 'e');

	setClippedEndPosition(gaps3, 3); //"-el"				//setClippedEndPosition
    SEQAN_ASSERT(gaps3 == "el");
	SEQAN_ASSERT_EQ(clippedEndPosition(gaps3), 3u);		//clippedEndPosition
	SEQAN_ASSERT_EQ(endPosition(gaps3), 3u);				//endPosition
    SEQAN_ASSERT_EQ(length(gaps3), 2u);

	setBeginPosition(gaps3, 0);	//"el"					//setBeginPosition
    SEQAN_ASSERT(gaps3[0] == 'e');
    SEQAN_ASSERT_EQ(beginPosition(gaps3), 0u);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gaps3), 1u);
    SEQAN_ASSERT_EQ(endPosition(gaps3), 2u);
    SEQAN_ASSERT_EQ(clippedEndPosition(gaps3), 3u);


//____________________________________________________________________________

	clear(gaps3);										//clear
	SEQAN_ASSERT_EQ(length(gaps3), 0u);	
    SEQAN_ASSERT_EQ(clippedBeginPosition(gaps3), 0u);
    SEQAN_ASSERT_EQ(clippedEndPosition(gaps3), 0u);

	setClippedEndPosition(gaps3, length(source(gaps3))); //reactivate after clear
    SEQAN_ASSERT(gaps3 == "hello");

	insertGaps(gaps3, 2, 3);
	setBeginPosition(gaps3, 2);
	SEQAN_ASSERT(gaps3 == "he---llo"); //"--he---llo"
    SEQAN_ASSERT_EQ(beginPosition(gaps3), 2u);

	//toSourcePosition
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 0), 0u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 1), 0u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 2), 0u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 3), 1u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 4), 2u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 5), 2u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 6), 2u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 7), 2u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 8), 3u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 9), 4u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 10), 5u);
    SEQAN_ASSERT_EQ(toSourcePosition(gaps3, 11), 5u);

	//toViewPosition
    SEQAN_ASSERT_EQ(toViewPosition(gaps3, 0), 2u);
    SEQAN_ASSERT_EQ(toViewPosition(gaps3, 1), 3u);
    SEQAN_ASSERT_EQ(toViewPosition(gaps3, 2), 7u);
    SEQAN_ASSERT_EQ(toViewPosition(gaps3, 3), 8u);
    SEQAN_ASSERT_EQ(toViewPosition(gaps3, 4), 9u);
    SEQAN_ASSERT_EQ(toViewPosition(gaps3, 5), 10u);

//____________________________________________________________________________

	SEQAN_ASSERT(gaps3 == "he---llo"); //"--he---llo"
    SEQAN_ASSERT_EQ(beginPosition(gaps3), 2u);

	clearGaps(gaps3, 1, 5);
	SEQAN_ASSERT(gaps3 == "he--llo"); //"-he--llo"
    SEQAN_ASSERT_EQ(beginPosition(gaps3), 1u);

	clearGaps(gaps3, 1, 5);
	SEQAN_ASSERT(gaps3 == "hello"); //"-hello"
    SEQAN_ASSERT_EQ(beginPosition(gaps3), 1u);

	clearGaps(gaps3);
	SEQAN_ASSERT(gaps3 == "hello"); //"hello"
    SEQAN_ASSERT_EQ(beginPosition(gaps3), 0u);

//____________________________________________________________________________

	SEQAN_ASSERT(gaps3 == "hello"); //"hello"

	setClippedBeginPosition(gaps3, 1);
	setClippedEndPosition(gaps3, 3); //"el"
    SEQAN_ASSERT(gaps3 == "el");
	SEQAN_ASSERT(sourceSegment(gaps3) == "el"); //sourceSegment
    SEQAN_ASSERT(sourceSegment(c_gaps3) == "el");


//____________________________________________________________________________
// Comparison Functions

	SEQAN_ASSERT(gaps3 == "el"); //"hello"
    SEQAN_ASSERT(gaps3 != "ello");
    SEQAN_ASSERT(gaps3 <= "el");
    SEQAN_ASSERT(gaps3 < "ello");
    SEQAN_ASSERT(gaps3 > "a");
    SEQAN_ASSERT(gaps3 >= "el");


//____________________________________________________________________________
}

//////////////////////////////////////////////////////////////////////////////

// Iterator Functions
template <typename TSource, typename TSpec>
void TestGapsIterator()
{
//____________________________________________________________________________
	typedef Gaps<TSource, TSpec> TGaps;


	typedef typename Iterator<TGaps, Rooted>::Type TIterator;

	TSource src1 = "hello";
	TGaps gaps4(src1);

	TIterator it1 = begin(gaps4);						//begin
	TIterator const & c_it1 = it1; //const version
	SEQAN_ASSERT(*it1 == 'h');							//operator *
	SEQAN_ASSERT(source(it1) == begin(src1));			//source
    SEQAN_ASSERT(source(c_it1) == begin(src1)) ;
	SEQAN_ASSERT(atBegin(c_it1));						//atBegin

	++it1;												//operator ++
    SEQAN_ASSERT(*it1 == 'e'); ;

	--it1;												//operator --
    SEQAN_ASSERT(*it1 == 'h'); ;

	TIterator it2 = end(gaps4);							//end
	SEQAN_ASSERT(atEnd(it2));							//atEnd

	--it2;
	SEQAN_ASSERT(*it2 == 'o');
//____________________________________________________________________________


	TIterator it3;										//default ctor
	TIterator it4 = it1;								//copy ctor
	TIterator const & c_it4 = it4;	

    SEQAN_ASSERT(container(it4) == container(it1));
    SEQAN_ASSERT(*it4 == *it1);

	SEQAN_ASSERT(it4 == it1);							//operator ==
    SEQAN_ASSERT(it4 == c_it1);
    SEQAN_ASSERT(c_it4 == it1);
    SEQAN_ASSERT(c_it4 == c_it1);

	++it1;

	SEQAN_ASSERT(it4 != it1);							//operator !=
    SEQAN_ASSERT(it4 != c_it1);
    SEQAN_ASSERT(c_it4 != it1);
    SEQAN_ASSERT(c_it4 != c_it1);

	it4 = it2;											//operator =
    SEQAN_ASSERT(*it4 == *it2);

	TIterator it5(gaps4, 1, 1);							//special ctor
//____________________________________________________________________________



//____________________________________________________________________________
}

//////////////////////////////////////////////////////////////////////////////

// Manipulation of Gaps
template <typename TSource, typename TSpec>
void TestGapManipulation()
{
//____________________________________________________________________________
	typedef Gaps<TSource, TSpec> TGaps;

//inserting gaps
	TSource src1 = "hello";
	TGaps gaps5(src1);

	insertGaps(gaps5, 1, 2); //insert gap somewhere
	SEQAN_ASSERT(gaps5 != "hello");
	SEQAN_ASSERT(gaps5 == "h--ello");

	insertGaps(gaps5, 1, 1); //insert blank at the beginning of a gap
	SEQAN_ASSERT(gaps5 == "h---ello");

	insertGaps(gaps5, 4, 1); //insert blank at the end of a gap
	SEQAN_ASSERT(gaps5 == "h----ello");

	insertGap(gaps5, 8); //insert second gap
	SEQAN_ASSERT(gaps5 == "h----ell-o");

	insertGaps(gaps5, 0, 2); //insert at position 0
	SEQAN_ASSERT(gaps5 == "h----ell-o"); //note: leading and trailing gaps are not displayed
	SEQAN_ASSERT_EQ(beginPosition(gaps5), 2u);

	insertGaps(gaps5, 8, 2); //insert gap with beginPosition == 2
	SEQAN_ASSERT(gaps5 == "h----e--ll-o");
	SEQAN_ASSERT_EQ(length(gaps5), 12u);

	insertGaps(gaps5, 14, 1); //insert gap behind end. Nothing happens
	SEQAN_ASSERT(gaps5 == "h----e--ll-o");
	SEQAN_ASSERT_EQ(length(gaps5), 12u);

//counting gaps
	SEQAN_ASSERT(gaps5 == "h----e--ll-o"); // "--h----e--ll-o"
	SEQAN_ASSERT_EQ(beginPosition(gaps5), 2u);

	SEQAN_ASSERT_EQ(countGaps(gaps5, 0), 2u);
	SEQAN_ASSERT_EQ(countGaps(gaps5, 2), 0u);
	SEQAN_ASSERT_EQ(countGaps(gaps5, 3), 4u);
	SEQAN_ASSERT_EQ(countGaps(gaps5, 4), 3u);
	SEQAN_ASSERT_EQ(countGaps(gaps5, 6), 1u);
	SEQAN_ASSERT_EQ(countGaps(gaps5, 8), 2u);
	SEQAN_ASSERT_EQ(countGaps(gaps5, 20), 0u);

//removing gaps
	SEQAN_ASSERT(gaps5 == "h----e--ll-o"); // "--h----e--ll-o"
	SEQAN_ASSERT_EQ(beginPosition(gaps5), 2u);

	removeGap(gaps5, 4); //remove gap somewhere in gap area
	SEQAN_ASSERT(gaps5 == "h---e--ll-o");

	removeGaps(gaps5, 6, 1); //try to remove gap in non-gap area. Nothing happens
	SEQAN_ASSERT("h---e--ll-o" == gaps5);

	removeGaps(gaps5, 7, 2); //remove gap region completely
	SEQAN_ASSERT(gaps5 == "h---ell-o");

	removeGaps(gaps5, 4, 10); //remove rest of gap region
	SEQAN_ASSERT(gaps5 == "h-ell-o");
//____________________________________________________________________________

}
//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSpec>
void TestTrailingGaps()
{
	typedef Gaps<TSource, TSpec> TGaps;

	//inserting gaps
	TSource src1 = "hello";
	TGaps gaps6(src1);
	SEQAN_ASSERT(gaps6 == "hello");

	insertGaps(gaps6, 10, 2);//somewhere behind the last char
	SEQAN_ASSERT(gaps6 == "hello"); //nothing changes
	SEQAN_ASSERT_EQ(countGaps(gaps6, 5), 0u);
	SEQAN_ASSERT_EQ(countGaps(gaps6, 10), 0u);

	insertGaps(gaps6, 5, 3); //directly behind the last char: append intended trailing gap
	SEQAN_ASSERT(gaps6 == "hello"); //"hello---"
	SEQAN_ASSERT_EQ(countGaps(gaps6,5), 3u);
	SEQAN_ASSERT_EQ(countGaps(gaps6,6), 2u);

	insertGap(gaps6, 8);//expand trailing gaps
	SEQAN_ASSERT(gaps6 == "hello"); //"hello----"
	SEQAN_ASSERT_EQ(countGaps(gaps6,5), 4u);
	SEQAN_ASSERT_EQ(countGaps(gaps6,6), 3u);
	SEQAN_ASSERT_EQ(countGaps(gaps6,7), 2u);
	SEQAN_ASSERT_EQ(countGaps(gaps6,8), 1u);
	SEQAN_ASSERT_EQ(countGaps(gaps6,9), 0u);
	SEQAN_ASSERT_EQ(countGaps(gaps6,10), 0u);
	SEQAN_ASSERT_EQ(countGaps(gaps6,20), 0u);

	insertGaps(gaps6, 9, 2);//expand trailing gaps on last position
	SEQAN_ASSERT(gaps6 == "hello"); //"hello------"
	SEQAN_ASSERT_EQ(countGaps(gaps6,5), 6u);

//removing gaps
	SEQAN_ASSERT(gaps6 == "hello"); // "hello------"
	SEQAN_ASSERT_EQ(beginPosition(gaps6), 0u);

	removeGap(gaps6, 7); //remove gap somewhere in gap area
	SEQAN_ASSERT(gaps6 == "hello");//"hello-----"
	SEQAN_ASSERT_EQ(countGaps(gaps6, 5), 5u);
	SEQAN_ASSERT_EQ(countGaps(gaps6, 7), 3u);
	SEQAN_ASSERT_EQ(countGaps(gaps6, 9), 1u);
	SEQAN_ASSERT_EQ(countGaps(gaps6, 10), 0u);

	removeGaps(gaps6, 8, 15); //remove rest of gap region
	SEQAN_ASSERT(gaps6 == "hello"); //"hello---"
	SEQAN_ASSERT_EQ(countGaps(gaps6, 5), 3u);
	SEQAN_ASSERT_EQ(countGaps(gaps6, 7), 1u);
	SEQAN_ASSERT_EQ(countGaps(gaps6, 8), 0u);

	removeGaps(gaps6, 5, 3); //remove gap region completely
	SEQAN_ASSERT(gaps6 == "hello"); //"hello"
	SEQAN_ASSERT_EQ(countGaps(gaps6, 5), 0u);

	insertGaps(gaps6, 5, 6);
	SEQAN_ASSERT_EQ(countGaps(gaps6,5), 6u);

	assignSource(gaps6, "new");		//clear trailing gaps when assign new source
	SEQAN_ASSERT(gaps6 == "new");   //"new"
	SEQAN_ASSERT_EQ(countGaps(gaps6, 3), 0u);
	insertGaps(gaps6, 3, 10);
	SEQAN_ASSERT_EQ(countGaps(gaps6, 3), 10u);

	TSource src2 = "hello";
	setSource(gaps6, src2);			//clear trailing gaps when set new source
	SEQAN_ASSERT(gaps6 == "hello"); //"re"
	SEQAN_ASSERT_EQ(countGaps(gaps6, 5), 0u);

	insertGaps(gaps6, 5, 10);
	SEQAN_ASSERT_EQ(countGaps(gaps6, 5), 10u);

	TGaps gaps6a(gaps6); 		//copy all trailing gaps
	SEQAN_ASSERT(gaps6a == "hello");
	SEQAN_ASSERT_EQ(countGaps(gaps6a, 5), 10u);

	TGaps gaps6b = gaps6a;		//assign the trailing gaps
	insertGaps(gaps6b, 10, 5);
	SEQAN_ASSERT(gaps6b == "hello");
	SEQAN_ASSERT_EQ(countGaps(gaps6a, 5), 10u);
	SEQAN_ASSERT_EQ(countGaps(gaps6b, 5), 15u);

	clearGaps(gaps6a);				//remove all trailing gaps
	SEQAN_ASSERT_EQ(countGaps(gaps6a, 5), 0u);
	SEQAN_ASSERT_EQ(countGaps(gaps6b, 5), 15u);

	}

template <typename TSource, typename TGapSpec>
inline void
TestCountCharacters() {
	typedef Gaps<TSource, TGapSpec> TGap;

	TSource seq = "hello";
	TGap gap7(seq);

	SEQAN_ASSERT(gap7 == "hello");
	SEQAN_ASSERT_EQ(length(gap7), 5u);

	SEQAN_ASSERT_EQ(countCharacters(gap7, 0), 5u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 1), 4u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 0), 5u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 3), 2u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 5), 0u);

	insertGaps(gap7, 0, 2);
	SEQAN_ASSERT(gap7 == "hello"); //"--hello"
	SEQAN_ASSERT_EQ(length(gap7), 5u);

	SEQAN_ASSERT_EQ(countCharacters(gap7, 0), 0u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 2), 5u);

	insertGap(gap7, 4);
	SEQAN_ASSERT(gap7 == "he-llo"); //"--he-llo"
	SEQAN_ASSERT_EQ(length(gap7), 6u);

	SEQAN_ASSERT_EQ(countCharacters(gap7, 0), 0u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 2), 2u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 3), 1u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 4), 0u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 5), 3u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 6), 2u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 7), 1u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 8), 0u);

	insertGaps(gap7, 6, 3);
	SEQAN_ASSERT(gap7 == "he-l---lo"); //"--he-l---lo"
	SEQAN_ASSERT_EQ(length(gap7), 9u);

	SEQAN_ASSERT_EQ(countCharacters(gap7, 0), 0u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 2), 2u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 3), 1u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 4), 0u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 5), 1u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 6), 0u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 7), 0u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 8), 0u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 9), 2u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 10), 1u);

	insertGaps(gap7, 11, 5);
	SEQAN_ASSERT(gap7 == "he-l---lo"); //"--he-l---lo-----"
	SEQAN_ASSERT_EQ(length(gap7), 9u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 11), 0u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 12), 0u);
	SEQAN_ASSERT_EQ(countCharacters(gap7, 30), 0u);
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_align_gaps_base_char_string_array_gaps) {
    TestGapsBase<String<char>, ArrayGaps>();
}


// SEQAN_DEFINE_TEST(test_align_gaps_base_char_string_sumlist_gaps) {
//     TestGapsBase<String<char>, SumlistGaps>();
// }


SEQAN_DEFINE_TEST(test_align_gaps_test_gaps_iterator) {
    TestGapsIterator<String<char>, ArrayGaps>();
}


SEQAN_DEFINE_TEST(test_align_gaps_test_gap_manipulation_char_string_array_gaps) {
	TestGapManipulation<String<char>, ArrayGaps>();
}


// SEQAN_DEFINE_TEST(test_align_gaps_test_gap_manipulation_char_string_sumlist_gaps) {
//     TestGapManipulation<String<char>, SumlistGaps>(); 
// }

SEQAN_DEFINE_TEST(test_align_gaps_test_trailing_gaps_char_string_array_gaps) {
	TestTrailingGaps<String<char>, ArrayGaps>();
}

SEQAN_DEFINE_TEST(test_align_gaps_test_count_characters_char_string_array_gaps) {
	TestCountCharacters<String<char>, ArrayGaps >();
}


#endif
