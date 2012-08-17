cmake_minimum_required (VERSION 2.6)
project (core_tests_align)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_align test_align.cpp
                                     test_align_gaps.h
                                     test_align_gaps_iterator.h
                                     test_align_global_alignment.h
                                     test_align_global_alignment_banded.h
                                     test_align_global_alignment_specialized.h
                                     test_alignment_base.h
                                     test_alignment_dp_band.h
                                     test_alignment_dp_value.h
                                     test_alignment_dp_matrix.h
                                     test_alignment_dp_formula.h
                                     test_alignment_dp_manager.h
                                     test_alignment_dp_tracker.h
                                     test_alignment_dp_impl.h
                                     test_alignment_traceback_tracesegment.h
                                     test_alignment_traceback_adaptor.h
                                     test_alignment_traceback_impl.h)
