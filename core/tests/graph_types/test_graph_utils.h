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

#ifndef SEQAN_HEADER_TEST_GRAPH_UTILS_H
#define SEQAN_HEADER_TEST_GRAPH_UTILS_H

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

template <typename TGraph>
void  Test_GraphDrawing_Tmp() {
	// Create a dummy graph
	TGraph g;
	addVertex(g);addVertex(g);addEdge(g,0,1);
	addVertex(g);addEdge(g,0,2);addVertex(g); 
	addVertex(g);addVertex(g);addEdge(g,3,4);

	// Dot Drawing
	write(std::cout,g,DotDrawing());
}

//////////////////////////////////////////////////////////////////////////////

void  Test_GraphDrawing() {
	typedef Position<String<char> >::Type TPosition;

	Test_GraphDrawing_Tmp<Graph<Directed<> > >();
	Test_GraphDrawing_Tmp<Graph<Undirected<> > >();
	Test_GraphDrawing_Tmp<Graph<Tree<> > >();

	// Automat
	Graph<Automaton<Dna> > automat;
	createRoot(automat);addVertex(automat);addEdge(automat,0,1,'a');
	addVertex(automat);addEdge(automat,0,2,'g');
	// Dot Drawing
	write(std::cout,automat,DotDrawing());

	// Trie
	Graph<Automaton<char> > trie;
	String<String<TPosition> > pos;
	String<String<char> > keywords;
	appendValue(keywords, String<char>("announce"));
	appendValue(keywords, String<char>("annual"));
	appendValue(keywords, String<char>("annually"));
	createTrie(trie,pos,keywords);
	// Dot Drawing
	String<String<char> > nodeMap;
	_createTrieNodeAttributes(trie, pos, nodeMap);
	String<String<char> > edgeMap;
	_createEdgeAttributes(trie,edgeMap);
	write(std::cout,trie,nodeMap, edgeMap, DotDrawing());

	// WordGraph
	typedef Graph<Automaton<Dna, String<Dna>, WordGraph<> > > TWordGraph;
	TWordGraph wordGr;
	addVertex(wordGr);addVertex(wordGr);addVertex(wordGr);addVertex(wordGr);
	assignRoot(wordGr,3);root(wordGr) = 2;addEdge(wordGr,0,3,"ag");
	addVertex(wordGr);addVertex(wordGr);addEdge(wordGr,0,5,"g");
	addVertex(wordGr);addVertex(wordGr);addEdge(wordGr,3,1,"aggg");
	addEdge(wordGr,3,4,"gg");addEdge(wordGr,5,2,"aggg");addEdge(wordGr,5,7,"g");
	addEdge(wordGr,7,6,"g");assignRoot(wordGr,0);
	write(std::cout,wordGr,DotDrawing());
}

	
//////////////////////////////////////////////////////////////////////////////


void Test_GraphUtils() {
	Test_GraphDrawing();
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_graph_utils)
{
	Test_GraphUtils();
}

#endif

