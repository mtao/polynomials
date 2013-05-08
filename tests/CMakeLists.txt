ENABLE_TESTING()

INCLUDE_DIRECTORIES(. ../include)
ADD_EXECUTABLE(initTest init_test.cpp)
ADD_EXECUTABLE(compactifyTest compactify_test.cpp)
ADD_EXECUTABLE(calculusTest calculus_test.cpp)


ADD_TEST(initTest ${CMAKE_CURRENT_BINARY_DIR}/initTest)
ADD_TEST(compactifyTest ${CMAKE_CURRENT_BINARY_DIR}/compactifyTest)
ADD_TEST(calculusTest ${CMAKE_CURRENT_BINARY_DIR}/calculusTest)



add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS initTest compactifyTest calculusTest)
