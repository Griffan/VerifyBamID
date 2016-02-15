cmake_minimum_required(VERSION 2.8)
project(ContaminationFinder)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp ContaminationEstimator.cpp MathGenMin.cpp MathGold.cpp Random.cpp SimplePileupViewer.cpp params.cpp)

add_subdirectory(statgen)
add_subdirectory(samtools)

include_directories(statgen samtools/htslib samtools)

LINK_DIRECTORIES(/Users/fanzhang/BioLib/lib/)

add_executable(ContaminationFinder ${SOURCE_FILES})

target_link_libraries(ContaminationFinder statgen samtools hts z)
#target_include_directories(ContaminationFinder PRIVATE "general")
#target_include_directories(ContaminationFinder PRIVATE "samtools")
#target_include_directories(ContaminationFinder PRIVATE "/Users/fanzhang/BioLib/lib/")