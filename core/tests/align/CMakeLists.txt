cmake_minimum_required (VERSION 2.6)
project (core_tests_align)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_align test_align.cpp
                                     test_align_align.h
                                     test_align_gaps.h
                                     test_align_matrix.h
                                     test_align_myers.h
                                     test_local_align.h)
