cmake_minimum_required(VERSION 3.14)
project(MoE)
find_package(Armadillo CONFIG REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 14)
add_executable(MoE main.cpp Expert.h Expert.cpp ExpertModel.h ExpertModel.cpp NormalModel.h NormalModel.cpp GLMModel.h GLMModel.cpp NormalFamily.h NormalFamily.cpp BinomialFamily.h BinomialFamily.cpp PoissonFamily.h PoissonFamily.cpp GammaFamily.h GammaFamily.cpp)
target_link_libraries(MoE ${ARMADILLO_LIBRARIES})
