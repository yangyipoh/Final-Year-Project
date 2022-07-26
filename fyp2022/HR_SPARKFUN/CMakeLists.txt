cmake_minimum_required(VERSION 3.1)
project (TBCode)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

find_package(pigpio REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(Threads REQUIRED)
# find_package(FFTW REQUIRED)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/wiringPi/wiringPi)
# include_directories(/home/pi/Arma/armadillo-8.600.0/include/armadillo_bits)
include_directories(/usr/local/lib)


set( NAME_SRC
    src/main.cpp 
    src/ppg.cpp
    src/MAX30105.cpp
    src/MPU.cpp
)

set(REQUIRED_LIBS
wiringPi
pthread
pigpio
lapack
)

add_executable(TBCode ${NAME_SRC})
target_link_libraries( TBCode ${REQUIRED_LIBS})



