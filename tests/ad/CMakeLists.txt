include(GoogleTest)

# Tests using Google Test
set(GTEST_TARGETS
    test_gemm
    test_greenstrain
    test_hadamard
    test_isotropic
    test_matdet
    test_matsum
    test_matinv
    test_mattrace
    test_matvecmult
    test_quaternion
    test_symeigs
    test_matmulttrace
    test_symrk
    test_symsum
    test_veccross
    test_vecnorm
    test_vecouter
    test_vecsum
    test_a2dmat
    test_a2dmatinv
    test_a2dmatdet
    test_ad_expressions
    test_adscalar
)

foreach(target IN LISTS GTEST_TARGETS)
    add_executable(${target} ${target}.cpp)

    target_include_directories(${target} PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/tests
    )

    target_link_libraries(${target} PRIVATE gtest_main)

    gtest_discover_tests(${target})
endforeach()

# # Add targets
# add_executable(test_ad_expressions test_ad_expressions.cpp)
# add_executable(test_a2dmat test_a2dmat.cpp)
# add_executable(test_a2dmatinv test_a2dmatinv.cpp)
# add_executable(test_a2dmatdet test_a2dmatdet.cpp)
# add_executable(test_adscalar test_adscalar.cpp)

# target_compile_options(test_ad_expressions PRIVATE -fsanitize=address)
# target_link_options(test_ad_expressions PRIVATE -fsanitize=address)

# # include A2D and test headers
# target_include_directories(test_ad_expressions PRIVATE
#     ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)
# target_include_directories(test_a2dmat PRIVATE
#     ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)
# target_include_directories(test_a2dmatinv PRIVATE
#     ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)
# target_include_directories(test_a2dmatdet PRIVATE
#     ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)
# target_include_directories(test_adscalar PRIVATE
#     ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)

# # For tests implmented using gtest, link them to gtest
# target_link_libraries(test_a2dmat PRIVATE gtest_main)
# target_link_libraries(test_a2dmatinv PRIVATE gtest_main)
# target_link_libraries(test_a2dmatdet PRIVATE gtest_main)

# include(GoogleTest)
# gtest_discover_tests(test_a2dmat)
# gtest_discover_tests(test_a2dmatinv)
# gtest_discover_tests(test_a2dmatdet)

# # Add non-gtest tests manually so that ctest could recognize it's a test
# add_test(NAME test_ad_expressions COMMAND test_ad_expressions)
# add_test(NAME test_adscalar COMMAND test_adscalar)

add_subdirectory(core)
