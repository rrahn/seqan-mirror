cmake_minimum_required (VERSION 2.6)
project (extras_apps_read_analyzer)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(read_analyzer read_analyzer.cpp
                                   analyze_454.h
                                   analyze_illumina.h
                                   histogram.h
                                   read_analyzer.h)

