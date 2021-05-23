cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

project(training_data)

include_directories(/usr/local/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../bin)

find_library(GILSP_LIB gilsp HINTS /usr/local/lib/)

############################################################
# run_gillespie
############################################################

add_executable(run_gillespie src/run_gillespie.cpp src/run_job_gillespie.hpp src/params.hpp src/reactions.hpp)

target_link_libraries(run_gillespie PUBLIC ${GILSP_LIB})

############################################################
# run_tau_leaping
############################################################

add_executable(run_tau_leaping src/run_tau_leaping.cpp src/run_job_tau_leaping.hpp src/params.hpp src/reactions.hpp)

target_link_libraries(run_tau_leaping PUBLIC ${GILSP_LIB})

############################################################
# extra
############################################################

# If want to include install target
# install(TARGETS bmla_layer_1 RUNTIME DESTINATION bin)