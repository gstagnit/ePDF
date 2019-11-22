# build tests
add_executable(Evolution Evolution.cc)
target_link_libraries(Evolution ePDF)
add_test(Evolution Evolution)
