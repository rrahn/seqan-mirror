cmake_minimum_required (VERSION 2.6)
project (sandbox_krakau_apps_bs_mason)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(bs_mason bs_mason.cpp bs_mason.h)
    # simulate_454.h
    #		            simulate_454_base_calling.h
    #		            simulate_illumina.h
    #		            simulate_illumina_bs.h
    #		            simulate_sanger.h
    #		            store_config.h
    #		            util.h )

# Add app tests if Python interpreter could be found.
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_bs_mason COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)
