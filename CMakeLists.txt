cmake_minimum_required(VERSION 3.10)
project(matrix)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I.")

# -----------------------------

FILE(GLOB MPI_SOURCES src/*.c src/mpi/*.c)
FILE(GLOB MPI_HEADERS include/libmatrix/*.h)

add_library(mpi_matrix STATIC 
    ${MPI_SOURCES}
    ${MPI_HEADERS}
)

target_compile_definitions(mpi_matrix PRIVATE USE_MPI)
target_include_directories(mpi_matrix PRIVATE include)

# -----------------------------

FILE(GLOB OMP_SOURCES src/*.c src/omp/*.c)
FILE(GLOB OMP_HEADERS include/libmatrix/*.h)

add_library(omp_matrix STATIC 
    ${OMP_SOURCES}
    ${OMP_HEADERS}
)

find_package(OpenMP REQUIRED)
target_include_directories(omp_matrix PRIVATE include)

# -----------------------------

FILE(GLOB BARE_SOURCES src/*.c src/bare/*.c)
FILE(GLOB BARE_HEADERS include/libmatrix/*.h)

add_library(bare_matrix STATIC 
    ${BARE_SOURCES} 
    ${BARE_HEADERS}
)

target_include_directories(bare_matrix PRIVATE include)

# -----------------------------

add_executable(bench_mpi bin/bench.c)
target_link_libraries(bench_mpi PRIVATE mpi_matrix)
target_include_directories(bench_mpi PRIVATE include)

add_executable(bench_omp bin/bench.c)
target_link_libraries(bench_omp PRIVATE omp_matrix)
target_include_directories(bench_omp PRIVATE include)

add_executable(bench_bare bin/bench.c)
target_link_libraries(bench_bare PRIVATE bare_matrix)
target_include_directories(bench_bare PRIVATE include)

add_executable(demo bin/demo.c)
target_link_libraries(demo PRIVATE mpi_matrix)
target_link_libraries(demo PRIVATE omp_matrix)
target_link_libraries(demo PRIVATE bare_matrix)
target_include_directories(demo PRIVATE include)

add_executable(tests bin/tests.c)
target_link_libraries(tests PRIVATE mpi_matrix)
target_link_libraries(tests PRIVATE omp_matrix)
target_link_libraries(tests PRIVATE bare_matrix)
target_include_directories(tests PRIVATE include)
