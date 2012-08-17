// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_H_

// TODO(holtgrew): Usage of gapped value in align module is not consistent, need proxies in many places, reference not cleanly implemented everywhere yet.
// TODO(holtgrew): The Author: tag at the top has to be corrected in the headers of this module.
// TODO(holtgrew): Anchor Gaps must be integrated completely.
// TODO(holtgrew): Local alignments & Fragments don't work nicely together at the moment, multiLocalAlignments() needs an equivalent in the new align module.
// TODO(holtgrew): Align<>, AlignCol<> need some love and documentation.
// TODO(holtgrew): Gaps need better documentation.
// TODO(holtgrew): graph_align should go away.
// TODO(holtgrew): refinement should become graph_align and also host Graph<Alignment<>>
// TODO(holtgrew): graph_msa should become align_msa, or so, see whitepaper.
// TODO(holtgrew): The documentation and Tutorial need heavy updates, clipping alignments howto can go away.

// ============================================================================
// Prerequisites
// ============================================================================

#include <seqan/basic.h>
#include <seqan/modifier.h>  // ModifiedAlphabet<>.
#include <seqan/refinement.h>  // TODO(holtgrew): This module should become graph_align and we should not depend on it.

// TODO(holtgrew): Why not use priority queue from STL?
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>

// ============================================================================
// Support
// ============================================================================

#include <seqan/align/gapped_value_type.h>
#include <seqan/align/align_config.h>
#include <seqan/align/align_traceback.h>
#include <seqan/align/matrix_base.h>

// ============================================================================
// Gaps & Gaps Iterator Data Structures
// ============================================================================

#include <seqan/align/gaps_base.h>
#include <seqan/align/gaps_iterator_base.h>

#include <seqan/align/gaps_array.h>
#include <seqan/align/gaps_iterator_array.h>

#include <seqan/align/gap_anchor.h>
#include <seqan/align/gaps_anchor.h>
#include <seqan/align/gaps_iterator_anchor.h>

// ============================================================================
// Alignment Data Structures and Columns
// ============================================================================

#include <seqan/align/align_metafunctions.h>
#include <seqan/align/align_cols.h>
#include <seqan/align/align_base.h>

// ============================================================================
// Alignment Algorithm Implementations.
// ============================================================================

// Defines all general tags and structures that are used for the alignment
// algorithms.
#include <seqan/align/alignment_base.h>

// The standard dp algorithms are unified in a single implementation of
// the alignment algorithm. We need different structures to distinguish
// between different properties of the standard dp algorithms.
// The following algorithms are consolidated:
// * Needleman-Wunsch and Banded Needleman-Wunsch,
// * Gotoh and Banded Gotoh,
// * Smith-Waterman and Banded Smith-Waterman,
// * Waterman-Eggert and Banded Waterman-Eggert TODO(rmaerker): Write me!

// Operations, structs and the implementation of the DP grid.
#include <seqan/align/alignment_dp_value.h>
#include <seqan/align/alignment_dp_band.h>
#include <seqan/align/alignment_dp_matrix.h>
#include <seqan/align/alignment_dp_formula.h>
#include <seqan/align/alignment_dp_tracker.h>
#include <seqan/align/alignment_dp_manager.h>
#include <seqan/align/alignment_dp_impl.h>
#include <seqan/align/alignment_dp.h>

// Operations and structures for computing the trace back and
// supporting different structures to store the alignment.
#include <seqan/align/alignment_traceback_tracesegment.h>
#include <seqan/align/alignment_traceback_impl.h>
#include <seqan/align/alignment_traceback_adaptor.h>
#include <seqan/align/alignment_traceback.h>

// The following algorithms could not be unified, so we list them
// separately.

// The Hirschberg algorithm computes an alignment including the trace back in
// linear space.
#include <seqan/align/global_alignment_hirschberg_impl.h>

// The implementations of Myers' bitvector algorithm for alignments can only
// compute alignment scores.  The combination of Hirschberg's and Myers'
// algorithm is limited in the same way.
#include <seqan/align/global_alignment_myers_impl.h>
#include <seqan/align/global_alignment_myers_hirschberg_impl.h>

// TODO(rmaerker): unify like other dp algorithms
// Implementations of the local alignment algorithms with declumping.  We also
// use them for the localAlignment() calls and return the best local alignment
// only.
#include <seqan/align/local_alignment_waterman_eggert_impl.h>
#include <seqan/align/local_alignment_banded_waterman_eggert_impl.h>

// ============================================================================
// Alignment Algorithm Interfaces
// ============================================================================

// Documentation for alignment algorithms.
#include <seqan/align/alignment_algorithm_interface.h>

// The front-end functions for global alignments.
#include <seqan/align/global_alignment_unbanded.h>
#include <seqan/align/global_alignment_banded.h>

// The front-end functions for local alignments.
#include <seqan/align/local_alignment_unbanded.h>
#include <seqan/align/local_alignment_banded.h>

// The front-end for enumeration of local alignments.
#include <seqan/align/local_alignment_enumeration.h>  // documentation
#include <seqan/align/local_alignment_enumeration_unbanded.h>
#include <seqan/align/local_alignment_enumeration_banded.h>

// The front-end functions for the more specialized alignment algorithms such as
// Hirschberg, Myers and Myers-Hirschberg.
#include <seqan/align/global_alignment_specialized.h>

// ============================================================================
// Operations On Alignments
// ============================================================================

#include <seqan/align/alignment_operations.h>

#endif  // SEQAN_CORE_INCLUDE_SEQAN_ALIGN_H_
