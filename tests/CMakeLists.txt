function(add_unittest TEST_NAME)
    set(TEST_EXE "unit${TEST_NAME}")
    # ${ARGN} stand for all other arguments (if any).
    add_executable(${TEST_EXE} "${ROOT_DIR}/tests/cubismup3d/${TEST_NAME}.cpp" ${ARGN})
    target_link_libraries(${TEST_EXE} cubismup3d::core)
    target_link_libraries(${TEST_EXE} cup_tests_lib)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_EXE})
endfunction()

add_library(cup_tests_lib
            "${ROOT_DIR}/tests/cubismup3d/Utils.cpp")
target_link_libraries(cup_tests_lib cubismup3d::core)
add_unittest(TestBoundaries)
add_unittest(TestBufferedLogger)
