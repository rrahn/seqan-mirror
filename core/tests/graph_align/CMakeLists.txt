cmake_minimum_required (VERSION 2.6)
project (core_tests_graph_align)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_graph_align test_graph_align.cpp
                                           test_graph_align.h)
