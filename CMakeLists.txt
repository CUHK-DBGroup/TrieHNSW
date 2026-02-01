cmake_minimum_required(VERSION 3.22)
project(LabelFilteredANNS)

set(CMAKE_CXX_STANDARD 17)

SET( CMAKE_CXX_FLAGS  "-O3 -march=native -lrt -DHAVE_CXX0X -fpic -w -fopenmp -ftree-vectorize -funroll-loops -mfma -fopenmp-simd -ftree-vectorizer-verbose=0 -mavx -mavx512f" )

option(ENABLE_SSE "Enable SSE support" ON)

add_executable(LabelFilteredANNS main.cpp
        Trie.hpp
        utils.hpp)
