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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_align_gaps.h"
#include "test_align_gaps_iterator.h"
// TODO(holtgrew): Test Align<>, AlignCols<>!
#include "test_align_global_alignment.h"
#include "test_align_global_alignment_banded.h"
#include "test_align_global_alignment_specialized.h"
#include "test_align_global_alignment_score.h"
#include "test_align_local_alignment.h"
#include "test_align_alignment_operations.h"

#include "test_alignment_base.h"
#include "test_alignment_dp_band.h"
#include "test_alignment_dp_value.h"
#include "test_alignment_dp_matrix.h"
#include "test_alignment_dp_formula.h"
#include "test_alignment_dp_manager.h"
#include "test_alignment_dp_tracker.h"
#include "test_alignment_dp_impl.h"
#include "test_alignment_traceback_tracesegment.h"
#include "test_alignment_traceback_adaptor.h"
#include "test_alignment_traceback_impl.h"

SEQAN_BEGIN_TESTSUITE(test_align)
{
    // -----------------------------------------------------------------------
    // Test Gaps Data Structures.
    // -----------------------------------------------------------------------

//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_metafunctions);
//    SEQAN_CALL_TEST(test_align_gaps_array_constructor_and_source);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_set_source);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_assign_source);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_gap_operations_gaps_center);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_gap_operations_gaps_leading);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_gap_operations_gaps_trailing);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_sequence_interface_ungapped);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_sequence_interface_gaps_center);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_sequence_interface_gaps_leading);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_sequence_interface_gaps_trailing);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_iterator_interface_begin);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_iterator_interface_end);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_iterator_interface_iter);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_source_view_position_ungapped);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_source_view_position_gaps_center);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_source_view_position_gaps_leading);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_source_view_position_gaps_trailing);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_clipping_ungapped);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_clipping_gaps_center);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_clipping_gaps_leading);
//    SEQAN_CALL_TEST(test_align_gaps_array_gaps_clipping_gaps_trailing);
//
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_metafunctions);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_constructor_and_source);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_set_source);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_assign_source);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_gap_operations_gaps_center);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_gap_operations_gaps_leading);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_gap_operations_gaps_trailing);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_sequence_interface_ungapped);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_sequence_interface_gaps_center);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_sequence_interface_gaps_leading);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_sequence_interface_gaps_trailing);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_iterator_interface_begin);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_iterator_interface_end);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_iterator_interface_iter);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_source_view_position_ungapped);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_source_view_position_gaps_center);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_source_view_position_gaps_leading);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_source_view_position_gaps_trailing);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_clipping_ungapped);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_clipping_gaps_center);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_clipping_gaps_leading);
//    SEQAN_CALL_TEST(test_align_gaps_anchor_gaps_clipping_gaps_trailing);
//
//    // -----------------------------------------------------------------------
//    // Test Gaps Iterators.
//    // -----------------------------------------------------------------------
//
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_metafunctions);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_trivial_iterator_array_functions);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_rooted_random_access_iterator_array_functions);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_movement);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_relations);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_pointer_arithmetic);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_forward_iteration);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_reverse_iteration);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_count_gaps_count_characters_is_gap);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_clipped_count_gaps_count_characters_is_gap);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_gap_operations_center);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_gap_operations_leading);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_array_gap_operations_trailing);
//
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_metafunctions);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_trivial_iterator_anchor_functions);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_rooted_random_access_iterator_anchor_functions);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_movement);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_relations);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_pointer_arithmetic);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_forward_iteration);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_reverse_iteration);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_count_gaps_count_characters_is_gap);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_clipped_count_gaps_count_characters_is_gap);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_gap_operations_center);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_gap_operations_leading);
//    SEQAN_CALL_TEST(test_align_gaps_iterator_anchor_gap_operations_trailing);
//
//
//    // -----------------------------------------------------------------------
//    // Test Alignment Structures and Methods
//    // -----------------------------------------------------------------------
//
//    // alignment_base.h
//    SEQAN_CALL_TEST(test_alignment_base_get_alignment_spec);
//    SEQAN_CALL_TEST(test_alignment_base_get_gap_spec);
//    SEQAN_CALL_TEST(test_alignment_base_is_traceback_on);
//    SEQAN_CALL_TEST(test_alignment_base_is_free_end_gap);
//    SEQAN_CALL_TEST(test_alignment_base_is_local);
//    SEQAN_CALL_TEST(test_alignment_base_is_global);
//
//    // alignment_dp_band.h
//    SEQAN_CALL_TEST(test_alignment_dp_band_constructor);
//    SEQAN_CALL_TEST(test_alignment_dp_band_get_horizontal_pos);
//    SEQAN_CALL_TEST(test_alignment_dp_band_get_vertical_pos);
//    SEQAN_CALL_TEST(test_alignment_dp_band_set_horizontal_pos);
//    SEQAN_CALL_TEST(test_alignment_dp_band_set_vertical_pos);
//    SEQAN_CALL_TEST(test_alignment_dp_band_get_band_size);
//    SEQAN_CALL_TEST(test_aligmenmt_dp_band_check_valid_band);
//    SEQAN_CALL_TEST(test_aligmenmt_dp_band_is_band_enabled);
//    SEQAN_CALL_TEST(test_aligmenmt_dp_band_is_wide_band);
//
//    // alignment_dp_value.h
//    SEQAN_CALL_TEST(test_alignment_dp_value_constructor);
//    SEQAN_CALL_TEST(test_alignment_dp_value_get_score);
//    SEQAN_CALL_TEST(test_alignment_dp_value_set_score);
//    SEQAN_CALL_TEST(test_alignment_dp_value_get_origin);
//    SEQAN_CALL_TEST(test_alignment_dp_value_set_origin);
//
//    // alingment_dp_matrix.h
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_constructor);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_value);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_reference);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_size);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_position);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_iterator);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_clear);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_length);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_begin);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_end);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_get_column_size);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_get_row_size);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_get_initial_column_size);
//    SEQAN_CALL_TEST(test_alignment_dp_matrix_get_initial_row_size);
//
//    // alignment_dp_formula.h
//    SEQAN_CALL_TEST(test_alignment_dp_formula_linear_diagonal);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_linear_vertical);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_linear_horizontal);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_linear_upper_band);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_linear_lower_band);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_linear_all);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_linear_none);
//
//    SEQAN_CALL_TEST(test_alignment_dp_formula_affine_diagonal);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_affine_vertical);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_affine_horizontal);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_affine_upper_band);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_affine_lower_band);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_affine_all);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_affine_none);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_global);
//    SEQAN_CALL_TEST(test_alignment_dp_formula_local);
//
//    // alignment_dp_manager
//    SEQAN_CALL_TEST(test_alignment_dp_manager_get_dp_direction);
//    SEQAN_CALL_TEST(test_alignment_dp_manager_is_to_track);
//    SEQAN_CALL_TEST(test_alignment_dp_manager_setup_column_manager);
//    SEQAN_CALL_TEST(test_alignment_dp_manager);
//    SEQAN_CALL_TEST(test_alignment_dp_manager_get_column_manager);
//    SEQAN_CALL_TEST(test_alignment_dp_manager_init);
//    SEQAN_CALL_TEST(test_alignment_dp_manager_update);
//    SEQAN_CALL_TEST(test_alignment_dp_manager_banded_chain);
//
//    // alignment_dp_tracker
//    SEQAN_CALL_TEST(test_alignment_dp_tracker_grid_coordinate);
//    SEQAN_CALL_TEST(test_alignment_dp_tracker_get_tracker_spec);
//    SEQAN_CALL_TEST(test_alignment_dp_tracker_go_begin_next_column);
//    SEQAN_CALL_TEST(test_alignment_dp_tracker_score_only);
//    SEQAN_CALL_TEST(test_alignment_dp_tracker_default);
////    SEQAN_CALL_TEST(test_alignment_dp_tracker_row_max);
////    SEQAN_CALL_TEST(test_alignment_dp_tracker_waterman_eggert);
//
//    // TODO (rmaerker): move to banded chain alignment algorithm
////    SEQAN_CALL_TEST(test_alignment_dp_tracker_banded_chain_init);
////    SEQAN_CALL_TEST(test_alignment_dp_tracker_banded_chain);
//
//    // alignment_dp_impl
//    SEQAN_CALL_TEST(test_alignment_dp_impl_begin_first_phase);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_begin_middle_phase);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_begin_last_phase);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_end_band);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_increment);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_trace);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_compute_cell);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_fill_column);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_is_valid_dp_settings);

    // Testing bands where upper diagonal crosses the horizontal sequence and the lower diagonal
    // crosses the vertical sequence.
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_noright_end_nottop_noleft_bottom_right);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_noright_end_nottop_noleft_bottom_noright);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_noright_end_nottop_noleft_nobottom_right);

    // Testing bands where upper diagonal crosses the horizontal sequence and the lower diagonal
    // crosses the horizontal sequence.
    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_noleft_nobottom_noright_end_nottop_noleft_bottom_right);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_noleft_nobottom_noright_end_nottop_noleft_bottom_noright);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_noleft_nobottom_noright_end_nottop_noleft_nobottom_right);
//
//    // Testing bands where upper diagonal crosses the horizontal sequence somewhere behind the
//    // end of the horizontal sequence and the lower diagonal crosses the horizontal sequence.
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_noleft_nobottom_right_end_top_noleft_bottom_right);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_noleft_nobottom_right_end_top_noleft_nobottom_right);

    // Testing bands where upper diagonal crosses the horizontal sequence somewhere behind the
    // end of the horizontal sequence and the lower diagonal crosses the vertical sequence.
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_right_end_top_noleft_bottom_right);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_left_nobottom_right_end_top_noleft_nobottom_right);
//
//    // Testing bands where upper diagonal crosses the vertical sequence and the lower diagonal
//    // crosses the vertical sequence.
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_notop_left_nobottom_norigth_end_notop_noleft_bottom_right);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_notop_left_nobottom_noright_end_notop_noleft_bottom_noright);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_notop_left_nobottom_noright_end_notop_noleft_nobottom_right);

    // Testing bands where upper diagonal crosses the vertical sequence and the lower diagonal
    // crosses the vertical sequence somewhere behind its end.
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_notop_left_nobottom_noright_end_notop_left_bottom_noright);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_notop_left_nobottom_noright_end_notop_left_bottom_right);
//
//    // Testing bands where upper diagonal crosses the horizontal sequence and the lower diagonal
//    // crosses the vertical sequence somewhere behind its end.
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_fullleft_end_bottom_noright);
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_top_fullleft_end_bottom_right);
//
//    // Testing bands where upper diagonal crosses the horizontal sequence somewhere behind its end
//    // and the lower diagonal crosses the vertical sequence somewhere behind its end.
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_fulltop_fullleft_end_fullbottom_fullright);
//
//    // Testing bands where upper diagonal crosses the horizontal sequence somewhere behind its end
//    // and the lower diagonal crosses the horizontal sequence somewhere behind its end.
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_begin_notop_noleft_nobottom_noright_end_notop_noleft_nobottom_noright);
//
//    // horizontal begin pos = vertical begin pos
//    SEQAN_CALL_TEST(test_alignment_dp_impl_band_one_size);

    // alignment_traceback_tracesegment.h
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_constructor);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_assignment);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_position);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_size);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_get_begin_horizontal);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_get_begin_vertical);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_get_end_horizontal);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_get_end_vertical);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_translate_trace_value);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_operator_stream);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_operator_equal);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_operator_unequal);
//    SEQAN_CALL_TEST(test_alignment_traceback_tracesegment_record_segment);
//
//    // alignment_traceback_adaptor.h
//    SEQAN_CALL_TEST(test_alignment_traceback_adaptor_adapt_file);
//    SEQAN_CALL_TEST(test_alignment_traceback_adaptor_adapt_gaps);
//    SEQAN_CALL_TEST(test_alignment_traceback_adaptor_adapt_fragments);
//    SEQAN_CALL_TEST(test_alignment_traceback_adaptor_adapt_alignment_graph);
//    SEQAN_CALL_TEST(test_alignment_traceback_adaptor_adapt_vertex_descriptor);
//
//    // alignment_traceback_impl.h
//
//    SEQAN_CALL_TEST(test_alignment_traceback_impl_follow_diagonal);
//    SEQAN_CALL_TEST(test_alignment_traceback_impl_follow_vertical);
//    SEQAN_CALL_TEST(test_alignment_traceback_impl_follow_horizontal);
//    SEQAN_CALL_TEST(test_alignment_traceback_impl_follow_trace);
//    SEQAN_CALL_TEST(test_alignment_traceback_impl_get_matrix_size);
//
//    // -----------------------------------------------------------------------
//    // Test Alignment Implementations
//    // -----------------------------------------------------------------------
//
//    // Unbanded Global Alignments
//
//    SEQAN_CALL_TEST(test_align_global_alignment_align_gaps_free_top_left_right_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_align_gaps_free_notop_left_noright_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_gaps_gaps_free_top_left_right_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_gaps_gaps_free_notop_left_noright_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_graph_gaps_free_top_left_right_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_graph_gaps_free_notop_left_noright_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_fragments_gaps_free_top_left_right_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_fragments_gaps_free_notop_left_noright_bottom_nw);
//
//    SEQAN_CALL_TEST(test_align_global_alignment_align_gaps_free_top_left_right_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_align_gaps_free_notop_left_noright_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_gaps_gaps_free_top_left_right_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_gaps_gaps_free_notop_left_noright_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_graph_gaps_free_top_left_right_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_graph_gaps_free_notop_left_noright_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_fragments_gaps_free_top_left_right_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_fragments_gaps_free_notop_left_noright_bottom_gotoh);
//
//    SEQAN_CALL_TEST(test_align_global_alignment_shorter_interfaces_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_shorter_interfaces_gotoh);
//
//    // -----------------------------------------------------------------------
//    // Test Banded Global Alignment Algorithms
//    // -----------------------------------------------------------------------
//
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_align_gaps_free_top_left_right_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_align_gaps_free_notop_left_noright_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_gaps_gaps_free_top_left_right_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_gaps_gaps_free_notop_left_noright_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_graph_gaps_free_top_left_right_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_graph_gaps_free_notop_left_noright_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_fragments_gaps_free_top_left_right_bottom_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_fragments_gaps_free_notop_left_noright_bottom_nw);
//
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_align_gaps_free_top_left_right_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_align_gaps_free_notop_left_noright_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_gaps_gaps_free_top_left_right_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_gaps_gaps_free_notop_left_noright_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_graph_gaps_free_top_left_right_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_graph_gaps_free_notop_left_noright_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_fragments_gaps_free_top_left_right_bottom_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_fragments_gaps_free_notop_left_noright_bottom_gotoh);
//
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_shorter_interfaces_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_banded_shorter_interfaces_gotoh);
//
//    // -----------------------------------------------------------------------
//    // Test Specialized Global Alignment Algorithms
//    // -----------------------------------------------------------------------
//
//    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_align);
//    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_gaps);
//    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_fragments);
//    SEQAN_CALL_TEST(test_align_global_alignment_hirschberg_graph);
//
//    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_align);
//    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_gaps);
//    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_fragments);
//    SEQAN_CALL_TEST(test_align_global_alignment_myers_hirschberg_graph);
//
//    // -----------------------------------------------------------------------
//    // Test Global Alignment Score Functions
//    // -----------------------------------------------------------------------
//
//    SEQAN_CALL_TEST(test_align_global_alignment_score_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_score_banded_nw);
//    SEQAN_CALL_TEST(test_align_global_alignment_score_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_score_banded_gotoh);
//    SEQAN_CALL_TEST(test_align_global_alignment_score_hirschberg);
//    SEQAN_CALL_TEST(test_align_global_alignment_score_myers);
//    SEQAN_CALL_TEST(test_align_global_alignment_score_myers_hirschberg);
//
//    // -----------------------------------------------------------------------
//    // Test Unbanded And Banded Local Alignment Algorithms
//    // -----------------------------------------------------------------------
//
//    SEQAN_CALL_TEST(test_align_local_alignment_align);
//    SEQAN_CALL_TEST(test_align_local_alignment_gaps);
//    SEQAN_CALL_TEST(test_align_local_alignment_graph);
//    SEQAN_CALL_TEST(test_align_local_alignment_fragment);
//
//    SEQAN_CALL_TEST(test_align_local_alignment_banded_align);
//    SEQAN_CALL_TEST(test_align_local_alignment_banded_gaps);
//    SEQAN_CALL_TEST(test_align_local_alignment_banded_graph);
//    SEQAN_CALL_TEST(test_align_local_alignment_banded_fragment);
//
//    // -----------------------------------------------------------------------
//    // Test Unbanded And Banded Local Alignment Enumeration Algorithms
//    // -----------------------------------------------------------------------
//
//    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_align);
//    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_gaps);
//    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_fragment);
//
//    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_banded_align);
//    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_banded_gaps);
//    SEQAN_CALL_TEST(test_align_local_alignment_enumeration_banded_fragment);

    // -----------------------------------------------------------------------
    // Test Operations On Align Objects
    // -----------------------------------------------------------------------

    SEQAN_CALL_TEST(test_align_integrate_align);
}
SEQAN_END_TESTSUITE
