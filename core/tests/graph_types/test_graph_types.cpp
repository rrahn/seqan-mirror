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

// External STL
#include <iostream>
#include <fstream>
#include <string>

// Seqan
#include <seqan/graph_types.h>

// Test files
#include "test_graph_basic.h"
#include "test_graph_types.h"
#include "test_graph_iterators.h"
#include "test_graph_properties.h"
#include "test_graph_derived.h"
#include "test_graph_utils.h"


SEQAN_BEGIN_TESTSUITE(test_graph_types) {
    // Call Tests.
    SEQAN_CALL_TEST(test_graph_basics);
    SEQAN_CALL_TEST(test_graph_types);	
    SEQAN_CALL_TEST(test_graph_iterators);
    SEQAN_CALL_TEST(test_graph_properties);
    SEQAN_CALL_TEST(test_graph_derived);
    SEQAN_CALL_TEST(test_graph_utils);
	
    // Verify Checkpoints.
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_interface.h");

	// basic
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_base.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_idmanager.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_edgestump.h");

	// types
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_impl_directed.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_impl_undirected.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_impl_automaton.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_impl_wordgraph.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_impl_tree.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_impl_fragment.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_impl_hmm.h");
	
	// iterators
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_iterator.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_iterator_vertex.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_iterator_outedge.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_iterator_adjacency.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_iterator_edge.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_iterator_bfs.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_iterator_dfs.h");
	
	// properties
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_property.h");
	
	// derived
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_impl_oracle.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_impl_trie.h");
	
	// utils	
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_drawing.h");
	SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_types/graph_utility_parsing.h");
}
SEQAN_END_TESTSUITE