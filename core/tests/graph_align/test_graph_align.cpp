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

#define TEST_PATH "projects/tests/graph_align/"

#include <seqan/basic.h>
#include <seqan/map.h>
#include <seqan/graph_align.h>
#include "test_graph_align.h"

SEQAN_BEGIN_TESTSUITE(test_graph_align) {
    // global alignments
    SEQAN_CALL_TEST(test_graph_align_needleman_wunsch);
    SEQAN_CALL_TEST(test_graph_align_gotoh);
    SEQAN_CALL_TEST(test_graph_align_hirschberg);
    SEQAN_CALL_TEST(test_graph_align_allAgainstAll);
    SEQAN_CALL_TEST(test_graph_align_gotohVsBandedGotoh);
    
    // local alignments
    SEQAN_CALL_TEST(test_graph_align_smith_waterman);
    SEQAN_CALL_TEST(test_graph_align_smith_waterman_clump);
    SEQAN_CALL_TEST(test_graph_align_banded_smith_waterman_clump);
    
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_align/graph_align_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_align/graph_align_config.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_align/graph_align_interface.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_align/graph_align_needleman_wunsch.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_align/graph_align_gotoh.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_align/graph_align_hirschberg.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_align/graph_align_smith_waterman.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/graph_align/graph_align_smith_waterman_clump.h");
}
SEQAN_END_TESTSUITE
