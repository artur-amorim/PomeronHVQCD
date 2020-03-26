# Specify the minimum CMAKE version required
cmake_minimum_required(VERSION 2.8.0)

set(CMAKE_CXX_FLAGS "-std=c++11")
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

# Fit HVQCD Global
add_executable(fitHVQCDGlobal.exe physics/fitHVQCDGlobal.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitHVQCDGlobal.exe armadillo)

# Fit HVQCD Without Scalars
add_executable(fitHVQCDWithoutScalars.exe physics/fitHVQCDWithoutScalars.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitHVQCDWithoutScalars.exe armadillo)

# Fit HVQCD Global sc and VgIR fixed
add_executable(fitHVQCDGlobal_sc_VgIR_fixed.exe physics/fitHVQCDGlobal_sc_VgIR_fixed.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitHVQCDGlobal_sc_VgIR_fixed.exe armadillo)

# Fit HVQCD Without Scalars sc and VgIR fixed
add_executable(fitHVQCDWithoutScalars_sc_VgIR_fixed.exe physics/fitHVQCDWithoutScalars_sc_VgIR_fixed.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitHVQCDWithoutScalars_sc_VgIR_fixed.exe armadillo)

# Fit HVQCD Without Glueball sc and VgIR fixed
add_executable(fitHVQCDWithoutGlueball_sc_VgIR_fixed.exe physics/fitHVQCDWithoutGlueball_sc_VgIR_fixed.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitHVQCDWithoutGlueball_sc_VgIR_fixed.exe armadillo)

# Fit HVQCD Without Glueball and Scalars sc and VgIR fixed
add_executable(fitHVQCDWithoutGlueballScalars_sc_VgIR_fixed.exe physics/fitHVQCDWithoutGlueballScalars_sc_VgIR_fixed.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitHVQCDWithoutGlueballScalars_sc_VgIR_fixed.exe armadillo)

# Fit HVQCD Without Glueball and Scalars
add_executable(fitHVQCDWithoutGlueballScalars.exe physics/fitHVQCDWithoutGlueballScalars.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitHVQCDWithoutGlueballScalars.exe armadillo)

# Fit HVQCD Without Glueball
add_executable(fitHVQCDWithoutGlueball.exe physics/fitHVQCDWithoutGlueball.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitHVQCDWithoutGlueball.exe armadillo)

# Fit the Singlet Axial Vector Mesons
add_executable(fitSingletAVM.exe physics/fitSingletAVM.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(fitSingletAVM.exe armadillo)

# Print HVQCD mass ratios
add_executable(HVQCDRatios.exe physics/HVQCDRatios.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(HVQCDRatios.exe armadillo)

# Save background fields, potentials and schrodinger potentials of the fluctuations
add_executable(savePotentials.exe physics/savePotentials.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(savePotentials.exe armadillo)

# Add test_schrodinger executable
add_executable(test_schrodinger.exe  tests/test_schrodinger.cpp ${SCHRODINGER})
target_link_libraries(test_schrodinger.exe armadillo)

# Add test_HVQCD executable
add_executable(test_HVQCD.exe  tests/test_HVQCD.cpp ${HVQCD} ${SCHRODINGER})
target_link_libraries(test_HVQCD.exe armadillo)

# Add test_EfficientNelderMead
add_executable(test_NelderMead.exe  tests/test_NelderMead.cpp)