# Specify the minimum CMAKE version required
cmake_minimum_required(VERSION 2.8.0)

set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_CXX_FLAGS "-std=c++11")
set(OS_INCLUDE_DIRS "/usr/include" "/usr/local/include")
set(OS_LIB_DIRS "/usr/lib" "/usr/local/lib")

ENABLE_LANGUAGE(Fortran)
set(CMKAE_Fortran_FLAGS "-std=legacy"
)
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

add_library(hvqcd_lib ${YANGMILLS} ${HVQCD} ${SCHRODINGER})

# Fit background routines
# fitYangMills 
add_executable(fitYangMills.exe physics/fitYangMills.cpp)
target_link_libraries(fitYangMills.exe armadillo hvqcd_lib)

# Fit HVQCD Global
add_executable(fitHVQCDGlobal.exe physics/fitHVQCDGlobal.cpp)
target_link_libraries(fitHVQCDGlobal.exe armadillo hvqcd_lib)

# Fit HVQCD Global a1 and a2 fixed
add_executable(fitHVQCDGlobal_a1_a2_fixed.exe physics/fitHVQCDGlobal_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDGlobal_a1_a2_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Without Scalars
add_executable(fitHVQCDWithoutScalars.exe physics/fitHVQCDWithoutScalars.cpp)
target_link_libraries(fitHVQCDWithoutScalars.exe armadillo hvqcd_lib)

# Fit HVQCD Without Scalars with a1 and a2 fixed
add_executable(fitHVQCDWithoutScalars_a1_a2_fixed.exe physics/fitHVQCDWithoutScalars_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDWithoutScalars_a1_a2_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Global sc and VgIR fixed
add_executable(fitHVQCDGlobal_sc_VgIR_fixed.exe physics/fitHVQCDGlobal_sc_VgIR_fixed.cpp)
target_link_libraries(fitHVQCDGlobal_sc_VgIR_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Global sc, VgIR, a1 and a2 fixed
add_executable(fitHVQCDGlobal_sc_VgIR_a1_a2_fixed.exe physics/fitHVQCDGlobal_sc_VgIR_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDGlobal_sc_VgIR_a1_a2_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Without Scalars sc and VgIR fixed
add_executable(fitHVQCDWithoutScalars_sc_VgIR_fixed.exe physics/fitHVQCDWithoutScalars_sc_VgIR_fixed.cpp)
target_link_libraries(fitHVQCDWithoutScalars_sc_VgIR_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Without Scalars sc, VgIR, a1 and a2 fixed
add_executable(fitHVQCDWithoutScalars_sc_VgIR_a1_a2_fixed.exe physics/fitHVQCDWithoutScalars_sc_VgIR_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDWithoutScalars_sc_VgIR_a1_a2_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Without Glueball sc and VgIR fixed
add_executable(fitHVQCDWithoutGlueball_sc_VgIR_fixed.exe physics/fitHVQCDWithoutGlueball_sc_VgIR_fixed.cpp)
target_link_libraries(fitHVQCDWithoutGlueball_sc_VgIR_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Without Glueball and Scalars sc and VgIR fixed
add_executable(fitHVQCDWithoutGlueballScalars_sc_VgIR_fixed.exe physics/fitHVQCDWithoutGlueballScalars_sc_VgIR_fixed.cpp)
target_link_libraries(fitHVQCDWithoutGlueballScalars_sc_VgIR_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Without Glueball and Scalars sc, VgIR, a1 and a2 fixed
add_executable(fitHVQCDWithoutGlueballScalars_sc_VgIR_a1_a2_fixed.exe physics/fitHVQCDWithoutGlueballScalars_sc_VgIR_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDWithoutGlueballScalars_sc_VgIR_a1_a2_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Without Glueball and Scalars
add_executable(fitHVQCDWithoutGlueballScalars.exe physics/fitHVQCDWithoutGlueballScalars.cpp)
target_link_libraries(fitHVQCDWithoutGlueballScalars.exe armadillo hvqcd_lib)

# Fit HVQCD Without Glueball and Scalars
add_executable(fitHVQCDWithoutGlueballScalars_a1_a2_fixed.exe physics/fitHVQCDWithoutGlueballScalars_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDWithoutGlueballScalars_a1_a2_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Without Glueball
add_executable(fitHVQCDWithoutGlueball.exe physics/fitHVQCDWithoutGlueball.cpp)
target_link_libraries(fitHVQCDWithoutGlueball.exe armadillo hvqcd_lib)

# Fit HVQCD Without Glueball a1 and a2 fixed
add_executable(fitHVQCDWithoutGlueball_a1_a2_fixed.exe physics/fitHVQCDWithoutGlueball_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDWithoutGlueball_a1_a2_fixed.exe armadillo hvqcd_lib)

# Fit the Singlet Axial Vector Mesons
add_executable(fitSingletAVM.exe physics/fitSingletAVM.cpp)
target_link_libraries(fitSingletAVM.exe armadillo hvqcd_lib)

# Print HVQCD mass ratios
add_executable(HVQCDRatios.exe physics/HVQCDRatios.cpp)
target_link_libraries(HVQCDRatios.exe armadillo hvqcd_lib)

# Save background fields, potentials and schrodinger potentials of the fluctuations
add_executable(savePotentials.exe physics/savePotentials.cpp)
target_link_libraries(savePotentials.exe armadillo hvqcd_lib)

# Add test_schrodinger executable
add_executable(test_schrodinger.exe  tests/test_schrodinger.cpp)
target_link_libraries(test_schrodinger.exe armadillo hvqcd_lib)

# Add test_HVQCD executable
add_executable(test_HVQCD.exe  tests/test_HVQCD.cpp)
target_link_libraries(test_HVQCD.exe armadillo hvqcd_lib)

# Add test_EfficientNelderMead
add_executable(test_NelderMead.exe  tests/test_NelderMead.cpp)