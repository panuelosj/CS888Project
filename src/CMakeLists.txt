include_directories (include)

add_executable(CS888Project main.cpp)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (CS888Project Eigen3::Eigen)

target_link_libraries(CS888Project libs)

add_library(MainLib main.cpp)
