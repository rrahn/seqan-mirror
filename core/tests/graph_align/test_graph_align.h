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

#ifndef SEQAN_HEADER_TEST_GRAPH_ALIGNMENT_H
#define SEQAN_HEADER_TEST_GRAPH_ALIGNMENT_H

#include <seqan/basic.h>
#include <seqan/map.h>
#include <seqan/graph_align.h>

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_graph_align_needleman_wunsch) {
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	
	TStringSet str;
	TString str0("annual");	assignValueById(str, str0);
	TString str1("annealing"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(0,-1,-1,0);
	int score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "annual");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "anneal");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 7)), "ing");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	int score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 3, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 3, score2);
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	TFragmentString matches;
	int score3 = globalAlignment(matches, stringSet(g), score_type, NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(length(matches), 1u);
	SEQAN_ASSERT_EQ(label(matches[0], stringSet(g), 0), "annual");
	SEQAN_ASSERT_EQ(label(matches[0], stringSet(g), 1), "anneal");
	SEQAN_ASSERT_EQ(score3, score);
	//std::cout << g << std::endl;

	str[0] = "annealing";
	str[1] = "annual";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "annual");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "anneal");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 7)), "ing");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 3, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 3, score2);
	//std::cout << g << std::endl;

	str[0] = "ThisisGarfieldthecat";
	str[1] = "Garfield";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "Thisis");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 9)), "Garfield");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 3)), "Garfield");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 17)), "thecat");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 4u);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 6, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score,score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 6, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 6, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,true>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 12, score2);
	//std::cout << g << std::endl;

	str[0] = "Garfield";
	str[1] = "ThisisGarfieldthecat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "Thisis");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 9)), "Garfield");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 3)), "Garfield");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 17)), "thecat");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 4u);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 6, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 6, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 6, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,true>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,true,false>(), NeedlemanWunsch() );
	SEQAN_ASSERT_EQ(score + 12, score2);
	//std::cout << g << std::endl;

	str[0] = "cat";
	str[1] = "ThisisGarfieldthecat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "ThisisGarfieldthe");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 18)), "cat");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "cat");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	//std::cout << g << std::endl;

	str[0] = "ThisisGarfieldthecat";
	str[1] = "cat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "ThisisGarfieldthe");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 18)), "cat");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "cat");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	//std::cout << g << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_graph_align_gotoh) {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	
	TStringSet str;
	TString str0("ttagt");	assignValueById(str, str0);
	TString str1("ttgt"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(1,-1,-1,-2);
	int score = globalAlignment(g, score_type, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 2)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 3)), "gt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 2)), "gt");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 5u);
	int score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	int score3 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_ASSERT_EQ(score2, score3);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_ASSERT_EQ(score2, score3);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_ASSERT_EQ(score2, score3);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_ASSERT_EQ(score2, score3);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_ASSERT_EQ(score2, score3);
	//std::cout << g << std::endl;

	str[0] = "ttgt";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 2)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 3)), "gt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 2)), "gt");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 5u);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_ASSERT_EQ(score2, score3);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_ASSERT_EQ(score2, score3);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_ASSERT_EQ(score2, score3);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_ASSERT_EQ(score2, score3);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_ASSERT_EQ(score2, score3);
	//std::cout << g << std::endl;

	str[0] = "tagt";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "tagt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 1)), "tagt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "t");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 2, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	//std::cout << g << std::endl;

	str[0] = "ttagt";
	str[1] = "tagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "tagt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 1)), "tagt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "t");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 2, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	//std::cout << g << std::endl;

	str[0] = "ttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 4)), "t");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 2, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 2, score2);
	//std::cout << g << std::endl;

	str[0] = "ttag";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 4)), "t");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 2, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 2, score2);
	//std::cout << g << std::endl;

	str[0] = "cttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 1)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "c");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 5)), "t");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "ttag");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 4u);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 2, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 2, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,true,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 4, score2);
	//std::cout << g << std::endl;

	str[0] = "ttag";
	str[1] = "cttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 1)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "c");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 5)), "t");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "ttag");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 4u);
	score2 = globalAlignment(g, score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 2, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 2, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,true,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 4, score2);
	//std::cout << g << std::endl;

	str[0] = "cttccagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "c");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 1)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 3)), "cc");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 5)), "ag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 7)), "t");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 2)), "ag");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 7u);
	//std::cout << g << std::endl;

	str[0] = "ttag";
	str[1] = "cttccagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "c");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 1)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 3)), "cc");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 5)), "ag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 7)), "t");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 2)), "ag");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 7u);
	//std::cout << g << std::endl;

	str[0] = "ttattaa";
	str[1] = "aaa";
	assignStringSet(g, str);
	Score<int> score_type2 = Score<int>(10,-1,-1,-2);
	score = globalAlignment(g, score_type2, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 2)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 3)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 5)), "aa");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 1)), "aa");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 6u);
	score2 = globalAlignment(stringSet(g), score_type2, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score + 3, score2);
	score2 = globalAlignment(stringSet(g), score_type2, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_ASSERT_EQ(score, score2);
	//std::cout << g << std::endl;

	str[0] = "aaa";
	str[1] = "ttattaa";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type2, Gotoh());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 2)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 3)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 5)), "aa");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 1)), "aa");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 6u);
	//std::cout << g << std::endl;


	// Note: Depending on the used recursion formula, the gotoh algorithms can differ !!!
	// Gotoh: Vertical gap and subsequent horizontal gap allowed
	// Gotoh3: Vertical gap and subsequent horizontal gap is not allowed
	//Score<double> score_scheme = Score<double>(5,-4,-0.5,-2);
	//str[0] = "tttggttt";
	//str[1] = "tttccttt";
	//assignStringSet(g, str);
	//score = globalAlignment(std::cout, str, score_scheme, Gotoh());
	//std::cout << std::endl;
	//score = globalAlignment(std::cout, str, score_scheme, Gotoh3());
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_graph_align_hirschberg) {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;

	TStringSet str;
	TString str0("ttagt");	assignValueById(str, str0);
	TString str1("ttgt"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(1,-1,-1,-2);
	int score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 2)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 3)), "gt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 2)), "gt");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 5u);
	int score2 = globalAlignment(stringSet(g), score_type, Hirschberg() );
	SEQAN_ASSERT_EQ(score, score2);
	//std::cout << g << std::endl;

	str[0] = "ttgt";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 2)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 3)), "gt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 2)), "gt");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 5u);
	//std::cout << g << std::endl;

	str[0] = "tagt";
	str[1] = "attagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "tagt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 2)), "tagt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "at");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	//std::cout << g << std::endl;

	str[0] = "attagt";
	str[1] = "tagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "tagt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 2)), "tagt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "at");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	//std::cout << g << std::endl;

	str[0] = "ttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 4)), "t");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	//std::cout << g << std::endl;

	str[0] = "ttag";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 4)), "t");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 3u);
	//std::cout << g << std::endl;

	str[0] = "cttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 1)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "c");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 5)), "t");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "ttag");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 4u);
	//std::cout << g << std::endl;

	str[0] = "ttag";
	str[1] = "cttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 1)), "ttag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "c");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 5)), "t");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "ttag");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	SEQAN_ASSERT_EQ(numVertices(g), 4u);
	//std::cout << g << std::endl;

	str[0] = "cttccagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "c");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 1)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 3)), "cc");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 5)), "ag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 7)), "t");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 2)), "ag");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 7u);
	//std::cout << g << std::endl;

	str[0] = "ttag";
	str[1] = "cttccagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "c");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 1)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 3)), "cc");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 5)), "ag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 7)), "t");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 2)), "ag");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 7u);
	//std::cout << g << std::endl;

	str[0] = "ttattaa";
	str[1] = "aaa";
	assignStringSet(g, str);
	Score<int> score_type2 = Score<int>(10,-1,-1,-2);
	score = globalAlignment(g, score_type2, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 2)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 3)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 5)), "aa");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 1)), "aa");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 6u);
	//std::cout << g << std::endl;

	str[0] = "aaa";
	str[1] = "ttattaa";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type2, Hirschberg());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 0)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 2)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 3)), "tt");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 5)), "aa");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 1)), "aa");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	SEQAN_ASSERT_EQ(numVertices(g), 6u);
	//std::cout << g << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, bool TBottom>
void TestGotohVSBandedGotoh_(AlignConfig<TTop, TLeft, TRight, TBottom> ac) {
	typedef unsigned int TSize;
	typedef int TScore;

	mtRandInit();
	for(TSize i = 0; i < 10; ++i) {
		typedef Dna5Q TAlphabet;
		typedef String<TAlphabet> TSequence;
		
		TSize lenN = mtRand() % 10 + 1;
		TSize lenM = mtRand() % 10 + 1;
		TSequence dna1;
		TSequence dna2;
		for(TSize i = 0; i<lenN; ++i) appendValue(dna1, TAlphabet(mtRand() % ValueSize<TAlphabet>::VALUE));
		for(TSize j = 0; j<lenM; ++j) appendValue(dna2, TAlphabet(mtRand() % ValueSize<TAlphabet>::VALUE));
		//dna1 = "AAAT";
		//dna2 = "CAA";
		
		TSize len1 = length(dna1);
		TSize len2 = length(dna2);
		typedef StringSet<TSequence, Dependent<> > TStringSet;
		typedef Graph<Alignment<TStringSet, void> > TGraph;
		TStringSet str;
		appendValue(str, dna1);
		appendValue(str, dna2);

		int matchScore =  mtRand() % 10;
		int misMatchScore =  -1 * (int) (mtRand() % 10);
		int gapScore =  -1 * (int) (mtRand() % 10);
		if (gapScore > misMatchScore) gapScore = misMatchScore;
		int gapOpenScore =  -1 * (int) (mtRand() % 10);
		if (gapOpenScore > gapScore) gapOpenScore = gapScore;
		//matchScore = 1;
		//misMatchScore = -1;
		//gapScore = -1;
		//gapOpenScore = -4;

		Score<int> score_type = Score<int>(matchScore,misMatchScore,gapScore,gapOpenScore);
		typedef String<Fragment<> > TFragmentString;
		TFragmentString matches;
		globalAlignment(matches, str, score_type, ac, Gotoh());
		int lowDiag = length(dna1);
		int highDiag = -1 * length(dna2);
		typedef typename Iterator<TFragmentString, Standard>::Type TFragIter;
		TFragIter itFrag = begin(matches, Standard());
		TFragIter itFragEnd = end(matches, Standard());
		for(;itFrag != itFragEnd; ++itFrag) {
			int newDiff = (int) itFrag->begin1 - (int) itFrag->begin2;
			if (newDiff > highDiag) highDiag = newDiff;
			if (newDiff < lowDiag) lowDiag = newDiff;
		}
		itFrag = begin(matches, Standard());
		if (itFrag != itFragEnd) {
			int pos1 = (int) itFrag->begin1 + (int) itFrag->len;
			int pos2 = (int) itFrag->begin2 + (int) itFrag->len;
			if (highDiag < ((int) len1 - pos2)) highDiag = ((int) len1 - pos2);
			if (lowDiag > (pos1 - (int) len2)) lowDiag = (pos1 - (int) len2);
			itFrag = itFragEnd;
			--itFrag;
			if (highDiag < (int) itFrag->begin1) highDiag = (int) itFrag->begin1;
			if (lowDiag > -1 * (int) itFrag->begin2) lowDiag = -1 * (int) itFrag->begin2;
		} else {
			lowDiag = -1 * length(dna2);
			highDiag = length(dna1);
		}
		TGraph g1(str);
		globalAlignment(g1, score_type, ac, Gotoh());
		int sc1 = _pairWiseSumOfPairsScore(g1, score_type, ac);
		TGraph g2(str);
		globalAlignment(g2, score_type, ac, lowDiag, highDiag, BandedGotoh());
		int sc2 = _pairWiseSumOfPairsScore(g2, score_type, ac);
		TGraph g3(str);
		globalAlignment(g3, score_type, ac, -1 * length(dna2), length(dna1), BandedGotoh());
		int sc3 = _pairWiseSumOfPairsScore(g3, score_type, ac);
		if ((sc1 != sc2) || (sc2 != sc3)) {
			std::cerr << "Randomized test failed:" << std::endl;
			std::cerr << "Seq1: " << dna1 << std::endl;
			std::cerr << "Seq2: " << dna2 << std::endl;
			std::cerr << "AlignConfig: " << TTop << ',' << TLeft << ',' << TRight << ',' << TBottom << std::endl;
			std::cerr << "Scores: (Matchscore: " << matchScore << ", Mismatchscore: " << misMatchScore << ", Gapscore: " << gapScore << ", GapOpenscore: " << gapOpenScore << ')' << std::endl;
			std::cerr << g1 << std::endl;
			std::cerr << "Score: " << sc1 << std::endl;
			std::cerr << g2 << std::endl;
			std::cerr << "Score: " << sc2 << std::endl;
			std::cerr << "Diagonals: " << lowDiag << ',' << highDiag << std::endl;
			std::cerr << g3 << std::endl;
			std::cerr << "Score: " << sc3 << std::endl;
			exit(-1);
		}
	}
}


SEQAN_DEFINE_TEST(test_graph_align_gotohVsBandedGotoh) {
	TestGotohVSBandedGotoh_(AlignConfig<false,false,false,false>() );
	TestGotohVSBandedGotoh_(AlignConfig<false,false,false,true>() );
	TestGotohVSBandedGotoh_(AlignConfig<false,false,true,false>() );
	TestGotohVSBandedGotoh_(AlignConfig<false,false,true,true>() );
	TestGotohVSBandedGotoh_(AlignConfig<false,true,false,false>() );
	TestGotohVSBandedGotoh_(AlignConfig<false,true,false,true>() );
	TestGotohVSBandedGotoh_(AlignConfig<false,true,true,false>() );
	TestGotohVSBandedGotoh_(AlignConfig<false,true,true,true>() );
	TestGotohVSBandedGotoh_(AlignConfig<true,false,false,false>() );
	TestGotohVSBandedGotoh_(AlignConfig<true,false,false,true>() );
	TestGotohVSBandedGotoh_(AlignConfig<true,false,true,false>() );
	TestGotohVSBandedGotoh_(AlignConfig<true,false,true,true>() );
	TestGotohVSBandedGotoh_(AlignConfig<true,true,false,false>() );
	TestGotohVSBandedGotoh_(AlignConfig<true,true,false,true>() );
	TestGotohVSBandedGotoh_(AlignConfig<true,true,true,false>() );
	TestGotohVSBandedGotoh_(AlignConfig<true,true,true,true>() );
}


//////////////////////////////////////////////////////////////////////////////////


template<bool TTop, bool TLeft, bool TRight, bool TBottom>
void AllAgainstAll__(AlignConfig<TTop, TLeft, TRight, TBottom> ac) {
	typedef unsigned int TSize;
	typedef int TScore;

	mtRandInit();
	for(TSize i = 0; i < 10; ++i) {
		typedef Dna5Q TAlphabet;
		typedef String<TAlphabet> TSequence;
		
		TSize lenN = mtRand() % 10 + 1;
		TSize lenM = mtRand() % 10 + 1;
		TSequence dna1;
		TSequence dna2;
		for(TSize i = 0; i<lenN; ++i) appendValue(dna1, TAlphabet(mtRand() % ValueSize<TAlphabet>::VALUE));
		for(TSize j = 0; j<lenM; ++j) appendValue(dna2, TAlphabet(mtRand() % ValueSize<TAlphabet>::VALUE));
		//dna1 = "AAGGTCA";
		//dna2 = "AT";
		
		TSize len1 = length(dna1);
		TSize len2 = length(dna2);
		typedef StringSet<TSequence, Dependent<> > TStringSet;
		typedef Graph<Alignment<TStringSet, void> > TGraph;
		TStringSet str;
		appendValue(str, dna1);
		appendValue(str, dna2);

		int matchScore =  mtRand() % 10;
		int misMatchScore =  -1 * (int) (mtRand() % 10);
		int gapScore =  -1 * (int) (mtRand() % 10);
		if (gapScore > misMatchScore) gapScore = misMatchScore;
		//matchScore = 0;
		//misMatchScore = -4;
		//gapScore = 0;

		Score<int> score_type = Score<int>(matchScore,misMatchScore,gapScore,gapScore);
		typedef String<Fragment<> > TFragmentString;
		TFragmentString matches;
		globalAlignment(matches, str, score_type, ac, NeedlemanWunsch());
		int lowDiag = length(dna1);
		int highDiag = -1 * length(dna2);
		typedef typename Iterator<TFragmentString, Standard>::Type TFragIter;
		TFragIter itFrag = begin(matches, Standard());
		TFragIter itFragEnd = end(matches, Standard());
		for(;itFrag != itFragEnd; ++itFrag) {
			int newDiff = (int) itFrag->begin1 - (int) itFrag->begin2;
			if (newDiff > highDiag) highDiag = newDiff;
			if (newDiff < lowDiag) lowDiag = newDiff;
		}
		itFrag = begin(matches, Standard());
		if (itFrag != itFragEnd) {
			int pos1 = (int) itFrag->begin1 + (int) itFrag->len;
			int pos2 = (int) itFrag->begin2 + (int) itFrag->len;
			if (highDiag < ((int) len1 - pos2)) highDiag = ((int) len1 - pos2);
			if (lowDiag > (pos1 - (int) len2)) lowDiag = (pos1 - (int) len2);
			itFrag = itFragEnd;
			--itFrag;
			if (highDiag < (int) itFrag->begin1) highDiag = (int) itFrag->begin1;
			if (lowDiag > -1 * (int) itFrag->begin2) lowDiag = -1 * (int) itFrag->begin2;
		} else {
			lowDiag = -1 * length(dna2);
			highDiag = length(dna1);
		}
		TGraph g1(str);
		globalAlignment(g1, score_type, ac, NeedlemanWunsch());
		int sc1 = _pairWiseSumOfPairsScore(g1, score_type, ac);
		TGraph g2(str);
		globalAlignment(g2, score_type, ac, lowDiag, highDiag, BandedNeedlemanWunsch());
		int sc2 = _pairWiseSumOfPairsScore(g2, score_type, ac);
		TGraph g3(str);
		globalAlignment(g3, score_type, ac, -1 * length(dna2), length(dna1), BandedNeedlemanWunsch());
		int sc3 = _pairWiseSumOfPairsScore(g3, score_type, ac);
		TGraph g4(str);
		globalAlignment(g4, score_type, ac, Gotoh());
		int sc4 = _pairWiseSumOfPairsScore(g4, score_type, ac);
		TGraph g5(str);
		globalAlignment(g5, score_type, ac, lowDiag, highDiag, BandedGotoh());
		int sc5 = _pairWiseSumOfPairsScore(g5, score_type, ac);
		TGraph g6(str);
		globalAlignment(g6, score_type, ac, -1 * length(dna2), length(dna1), BandedGotoh());
		int sc6 = _pairWiseSumOfPairsScore(g6, score_type, ac);
		if ((sc1 != sc2) || (sc2 != sc3) || (sc3 != sc4) || (sc4 != sc5) || (sc5 != sc6)) {
			std::cerr << "Randomized test failed:" << std::endl;
			std::cerr << "Seq1: " << dna1 << std::endl;
			std::cerr << "Seq2: " << dna2 << std::endl;
			std::cerr << "AlignConfig: " << TTop << ',' << TLeft << ',' << TRight << ',' << TBottom << std::endl;
			std::cerr << "Scores: (Matchscore: " << matchScore << ", Mismatchscore: " << misMatchScore << ", Gapscore: " << gapScore << ')' << std::endl;
			std::cerr << g1 << std::endl;
			std::cerr << "Score: " << sc1 << std::endl;
			std::cerr << g2 << std::endl;
			std::cerr << "Score: " << sc2 << std::endl;
			std::cerr << "Diagonals: " << lowDiag << ',' << highDiag << std::endl;
			std::cerr << g3 << std::endl;
			std::cerr << "Score: " << sc3 << std::endl;
			std::cerr << g4 << std::endl;
			std::cerr << "Score: " << sc4 << std::endl;
			std::cerr << g5 << std::endl;
			std::cerr << "Score: " << sc5 << std::endl;
			std::cerr << "Diagonals: " << lowDiag << ',' << highDiag << std::endl;
			std::cerr << g6 << std::endl;
			std::cerr << "Score: " << sc6 << std::endl;
			exit(-1);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_graph_align_allAgainstAll) {
	AllAgainstAll__(AlignConfig<false,false,false,false>() );
	AllAgainstAll__(AlignConfig<false,false,false,true>() );
	AllAgainstAll__(AlignConfig<false,false,true,false>() );
	AllAgainstAll__(AlignConfig<false,false,true,true>() );
	AllAgainstAll__(AlignConfig<false,true,false,false>() );
	AllAgainstAll__(AlignConfig<false,true,false,true>() );
	AllAgainstAll__(AlignConfig<false,true,true,false>() );
	AllAgainstAll__(AlignConfig<false,true,true,true>() );
	AllAgainstAll__(AlignConfig<true,false,false,false>() );
	AllAgainstAll__(AlignConfig<true,false,false,true>() );
	AllAgainstAll__(AlignConfig<true,false,true,false>() );
	AllAgainstAll__(AlignConfig<true,false,true,true>() );
	AllAgainstAll__(AlignConfig<true,true,false,false>() );
	AllAgainstAll__(AlignConfig<true,true,false,true>() );
	AllAgainstAll__(AlignConfig<true,true,true,false>() );
	AllAgainstAll__(AlignConfig<true,true,true,true>() );
}


//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_graph_align_smith_waterman) {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;

	TStringSet str;
	TString str0("gctctgcgaata"); assignValueById(str, str0);
	TString str1("cgttgagatact"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(2,-1,-2,-2);
	int score = localAlignment(g, score_type, SmithWaterman());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 4)), "tgcg");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 8)), "a");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 9)), "ata");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 3)), "tgag");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 7)), "ata");
	SEQAN_ASSERT_EQ(numEdges(g), 2u);
	int score2 = localAlignment(stringSet(g), score_type, SmithWaterman());
	SEQAN_ASSERT_EQ(score, score2);
	//std::cout << g << std::endl;
	
	str[0] = "TTGACACCCTCCCAATTGTA";
	str[1] = "ACCCCAGGCTTTACACAT";
	assignStringSet(g, str);
	score = localAlignment(g, score_type, SmithWaterman());
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 0, 0)), "TTGACAC");
	SEQAN_ASSERT_EQ(label(g, findVertex(g, 1, 9)), "TTTACAC");
	SEQAN_ASSERT_EQ(numEdges(g), 1u);
	score2 = localAlignment(stringSet(g), score_type, SmithWaterman());
	SEQAN_ASSERT_EQ(score, score2);
	//std::cout << g << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_graph_align_smith_waterman_clump) {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;

	TStringSet str;
	TString str0("gctctgcgaata"); assignValueById(str, str0);
	TString str1("cgttgagatact"); assignValueById(str, str1);
	Score<int> score_type = Score<int>(2,-1,-2,-2);
	String<Fragment<> > matches;
	String<int> scores;
	multiLocalAlignment(str, matches, scores, score_type, 4, SmithWatermanClump() );
	//std::cout << length(matches) << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_graph_align_banded_smith_waterman_clump) {
    typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;

	TStringSet str;
	TString str0("ggggcttaagcttgggg"); assignValueById(str, str0);
	TString str1("aaaacttagctctaaaa"); assignValueById(str, str1);
	Score<int> score_type = Score<int>(2,-1,-2,-2);
	String<Graph<Alignment<TStringSet> > > alignments;
	String<int> scores;
	multiLocalAlignment(str, alignments, scores, score_type, 5, -6, 6, BandedSmithWatermanClump() );

   SEQAN_ASSERT_EQ(length(alignments), 4u);
   SEQAN_ASSERT_EQ(length(scores), 4u);

   SEQAN_ASSERT_EQ(value(scores, 0), 12);
   Graph<Alignment<TStringSet> > align = value(alignments, 0);
   SEQAN_ASSERT_EQ(label(align, findVertex(align, 0, 0)), "GGGG");
   SEQAN_ASSERT_EQ(label(align, findVertex(align, 0, 4)), "CTT");
   SEQAN_ASSERT_EQ(label(align, findVertex(align, 0, 7)), "A");
   SEQAN_ASSERT_EQ(label(align, findVertex(align, 0, 8)), "AGCT");
   SEQAN_ASSERT_EQ(label(align, findVertex(align, 0, 12)), "TGGGG");
   SEQAN_ASSERT_EQ(label(align, findVertex(align, 1, 0)), "AAAA");
   SEQAN_ASSERT_EQ(label(align, findVertex(align, 1, 4)), "CTT");
   SEQAN_ASSERT_EQ(label(align, findVertex(align, 1, 7)), "AGCT");
   SEQAN_ASSERT_EQ(label(align, findVertex(align, 1, 11)), "CTAAAA");
}

#endif

