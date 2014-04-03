execute_process(COMMAND perf stat -o /dev/stdout w OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE PERF-FOUND)
if (${PERF-FOUND} EQUAL 0)
  message(STATUS "perf was found and accepted the -o option")
else()
  message(STATUS "perf was not found or did not accept the -o option")
endif()

function(add_tests dir)
  # getting a list of tests
  execute_process(
    COMMAND "find" "${dir}" "-name" "*:in"
    COMMAND "sort"
    COMMAND "cut" "-d" ":" "-f" "1"
    COMMAND "awk" "{printf(\"%s;\",$1)}" # to make it a CMake list
    OUTPUT_VARIABLE TESTS
  )

  # registering the tests to be run by CTest 
  foreach(TEST ${TESTS})
    string(STRIP ${TEST} TEST)
    get_filename_component(NAME "${TEST}" NAME)

    # parsing filenames
    string(REGEX MATCHALL "nx=(.+)_ny=(.+)_Cx=(.+)_Cy=(.+)_nt=(.+)_it=(.+)" NX ${NAME})
    set(NX ${CMAKE_MATCH_1}) # \__ domain size
    set(NY ${CMAKE_MATCH_2}) # /
    set(CX ${CMAKE_MATCH_3}) # \__ Courant numbers
    set(CY ${CMAKE_MATCH_4}) # /
    set(NT ${CMAKE_MATCH_5}) # --- number of timesteps
    set(IT ${CMAKE_MATCH_6}) # --- number of MPDATA iterations
    set(DEC 4) # ----------------- desired precision (abs(desired-actual) < 0.5 * 10**(-decimal))


    set(PFX "--format=time:")
    set(SFX "\;nx=\;${NX}\;ny=\;${NY}\;nt=\;${NT}\;it=\;${IT}\;peakrss=\;%M\;cputime=\;%U\;exit=\;%x\;")

    # Linux perf
    set(PERFTIME "perf;stat;-o;/dev/stdout;--event")
    # 
    set(PERFTIME "${PERFTIME};cache-misses,cache-references")
    # L1 data cache
    set(PERFTIME "${PERFTIME},L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores,L1-dcache-store-misses,L1-dcache-prefetches,L1-dcache-prefetch-misses")
    # L1 instruction cache
    set(PERFTIME "${PERFTIME},L1-icache-loads,L1-icache-load-misses,L1-icache-prefetches,L1-icache-prefetch-misses")
    # last level cache
    set(PERFTIME "${PERFTIME},LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,LLC-prefetches,LLC-prefetch-misses")
    # data translation lookaside buffer
    set(PERFTIME "${PERFTIME},dTLB-loads,dTLB-load-misses,dTLB-stores,dTLB-store-misses,dTLB-prefetches,dTLB-prefetch-misses")
    # instruction translation lookaside buffer
    set(PERFTIME "${PERFTIME},iTLB-loads,iTLB-load-misses")
    # cycles and instructions
    set(PERFTIME "${PERFTIME},cycles,instructions")
    # GNU time
    if (${PERF-FOUND} EQUAL 0)
      set(PERFTIME "${PERFTIME};time")
    else()
      set(PERFTIME "time")
    endif()

    set(REPEAT "0;1;2;") # repeat three times
    set(VERS "gcc;cxx;pyt;pyp;for") # repeat for each version

    foreach(V ${VERS})

      # setting CMD
      if (${V} STREQUAL "gcc")
        # C++ - GCC
        set(CMD ${CMAKE_BINARY_DIR}/cpp/test-cpp-gnucxx)
      elseif(${V} STREQUAL "cxx")
        # C++ - non-GNU compiler (e.g. Clang++)
        set(CMD ${CMAKE_BINARY_DIR}/cpp/test-cpp-nongnu)
      elseif(${V} STREQUAL "pyt")
	# Python
        set(CMD python;-OO;${CMAKE_BINARY_DIR}/pyt/test.py)
      elseif(${V} STREQUAL "pyp")
        # PyPy
        set(CMD pypy;${CMAKE_BINARY_DIR}/pyt/test.py)
      elseif(${V} STREQUAL "for")
        # Fortran
        set(CMD ${CMAKE_BINARY_DIR}/for/test-for)
      else()
        message(FATAL_ERROR "fix VERS! (${V})") 
      endif()       

      # running it N-times using CMD for NT=0 (i.e. data load test)
      foreach(N ${REPEAT})
	set(TIME_ARG "${PFX}${V}-load-${N}${SFX}")
	set(ARG ${NX};${NY};${CX};${CY};0;${IT};${TEST}:in;${TEST}:in;${DEC})
	add_test(NAME "${V}:${NAME}-load-${N}" COMMAND time ${TIME_ARG} ${CMD} ${ARG})
      endforeach()

      # running it N-times using CMD for NT=NT 
      foreach(N ${REPEAT})
	set(TIME_ARG "${PFX}${V}-calc-${N}${SFX}")
	set(ARG ${NX};${NY};${CX};${CY};${NT};${IT};${TEST}:in;${TEST}:out;${DEC})
	add_test(NAME "${V}:${NAME}-calc-${N}" COMMAND ${PERFTIME} ${TIME_ARG} ${CMD} ${ARG})
      endforeach()

    endforeach()
  endforeach()
endfunction()

add_subdirectory(sanity)
add_subdirectory(timing)
