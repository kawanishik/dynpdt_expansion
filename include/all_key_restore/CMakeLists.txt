set(CMAKE_CXX_STANDARD 17)

enable_testing()

add_executable(test_all_key_restore test_all_key_restore.cpp)
add_test(test_all_key_restore test_all_key_restore)