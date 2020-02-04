# Specify the minimum CMAKE version required
cmake_minimum_required(VERSION 2.8.0)

project(PomeronHVQCD VERSION 0.0.1)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")
set(OS_INCLUDE_DIRS "/usr/include" "/usr/local/include")
set(OS_LIB_DIRS "/usr/lib" "/usr/local/lib")

include_directories(${OS_INCLUDE_DIRS} "include")
link_directories(${OS_LIB_DIRS})

make_directory(execs)

set(EXECUTABLE_OUTPUT_PATH "execs")

# Variable SCHRODINGER to be used when schrodinger src code is needed
set(SCHRODINGER src/schrodinger/common.cpp src/schrodinger/solvspec.cpp src/schrodinger/chebspec.cpp src/schrodinger/numerov.cpp src/schrodinger/schrodinger.cpp)

# Variable YangMills
set(YANGMILLS src/background.cpp src/YangMills.cpp)

# Variable HVQCD
set(HVQCD src/background.cpp src/HolographicVQCD.cpp)

# Fit background routines
# fitYangMills 
add_executable(fitYangMills.exe physics/fitYangMills.cpp ${YANGMILLS} ${SCHRODINGER})
target_link_libraries(fitYangMills.exe armadillo)

# add_executable(fitBackground.exe physics/fitBackground.cpp ${HVQCD} ${SCHRODINGER})
# target_link_libraries(fitBackground.exe armadillo)

# Print HVQCD mass ratios
#add_executable(HVQCD_mass_ratios.exe physics/HVQCD_mass_ratios.cpp ${HVQCD} ${SCHRODINGER})
#target_link_libraries(HVQCD_mass_ratios.exe armadillo)

# Add test_schrodinger executable
add_executable(test_schrodinger.exe  tests/test_schrodinger.cpp ${SCHRODINGER})
target_link_libraries(test_schrodinger.exe armadillo)

# Add test_HVQCD executable
add_executable(test_HVQCD.exe  tests/test_HVQCD.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(test_HVQCD.exe armadillo)