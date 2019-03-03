function(add_unittest TEST_NAME)
    set(TEST_EXE "unit${TEST_NAME}")
    add_executable(${TEST_EXE} ${ARGN})
    target_link_libraries(${TEST_EXE} ${STATIC_LIB})  # cubismup3d
    target_link_libraries(${TEST_EXE} cup_tests_lib)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_EXE})
endfunction()

add_library(cup_tests_lib
            ${ROOT_FOLDER}/tests/Utils.cpp)
add_unittest(TestBoundaries ${ROOT_FOLDER}/tests/cubismup3d/TestBoundaries.cpp)
