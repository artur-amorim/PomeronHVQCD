# Specify the minimum CMAKE version required
cmake_minimum_required(VERSION 3.14.4)

project(PomeronHVQCD VERSION 0.0.1)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3")
set(OS_INCLUDE_DIRS "/usr/include" "/usr/local/include")
set(OS_LIB_DIRS "/usr/lib" "/usr/local/lib")

include_directories(${OS_INCLUDE_DIRS} "include")
link_directories(${OS_LIB_DIRS})

# Variable SCHRODINGER to be used when schrodinger src code is needed
set(SCHRODINGER src/schrodinger/common.cpp src/schrodinger/solvspec.cpp src/schrodinger/chebspec.cpp src/schrodinger/numerov.cpp src/schrodinger/schrodinger.cpp)

# Variable HVQCD
set(HVQCD src/background.cpp src/HolographicVQCD.cpp)

# Fit background routines
add_executable(fitBackground_kU1_not_fixed.exe physics/fitBackground_kU1_not_fixed.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitBackground_kU1_not_fixed.exe armadillo)

add_executable(fitSingletAxialVectorMesons.exe physics/fitSingletAxialVectorMesons.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitSingletAxialVectorMesons.exe armadillo)

add_executable(fittingBackground.exe physics/fittingBackground.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fittingBackground.exe armadillo)

# Print HVQCD mass ratios
add_executable(HVQCD_mass_ratios.exe physics/HVQCD_mass_ratios.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(HVQCD_mass_ratios.exe armadillo)

# Add test_schrodinger executable
add_executable(test_schrodinger.exe  tests/test_schrodinger.cpp ${SCHRODINGER})
target_link_libraries(test_schrodinger.exe armadillo)

# Add test_HVQCD executable
add_executable(test_HVQCD.exe  tests/test_HVQCD.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(test_HVQCD.exe armadillo)