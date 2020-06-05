# Specify the minimum CMAKE version required
cmake_minimum_required(VERSION 2.8.0)

set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_CXX_FLAGS "-std=c++11 -O3")
ENABLE_LANGUAGE(Fortran)
set(CMAKE_Fortran_FLAGS "-std=legacy")

if(UNIX AND NOT APPLE)
        FIND_PACKAGE(Boost 1.68 REQUIRED)
        include_directories(${Boost_INCLUDE_DIRS} "include")
endif()

if(APPLE)
        set(OS_INCLUDE_DIRS "/usr/include" "/usr/local/include")
        set(OS_LIB_DIRS "/usr/lib" "/usr/local/lib")
        include_directories(${OS_INCLUDE_DIRS} "include")
        link_directories(${OS_LIB_DIRS})
endif()

make_directory(execs)

set(EXECUTABLE_OUTPUT_PATH "execs")

# Variable SCHRODINGER to be used when schrodinger src code is needed
set(SCHRODINGER src/schrodinger/common.cpp src/schrodinger/solvspec.cpp src/schrodinger/chebspec.cpp src/schrodinger/numerov.cpp src/schrodinger/matrixNumerov.cpp src/schrodinger/schrodinger.cpp)
#add_library(schrodinger_lib ${SCHRODINGER})

# QUADPACK
set(QUADPACK src/quadpack/d1mach.f src/quadpack/dqags.f src/quadpack/dqagse.f src/quadpack/dqelg.f src/quadpack/dqk21.f src/quadpack/dqpsrt.f src/quadpack/fdump.f src/quadpack/i1mach.f src/quadpack/j4save.f src/quadpack/s88fmt.f src/quadpack/xerabt.f src/quadpack/xerctl.f src/quadpack/xerprt.f src/quadpack/xerror.f src/quadpack/xerrwv.f src/quadpack/xersav.f src/quadpack/xgetua.f)
#add_library(quadpack_lib ${QUADPACK})

# COLNEW
set(COLNEW src/Fortran/colnew.f src/Fortran/dgefa.f src/Fortran/dgesl.f)
#add_library(colnew_lib ${COLNEW})

# Variable YangMills
set(YANGMILLS src/background.cpp src/YangMills.cpp)

# Variable HVQCD
set(HVQCD src/background.cpp src/HolographicVQCD.cpp)

# U1NNMode
set(U1NNMODE src/U1NNMode.cpp)

# Spectra and Reggeon classes
set(SPECTRA src/Spectra.cpp src/Reggeon.cpp src/Kernel.cpp src/GluonKernel.cpp)

# ProcessObservable and its derived classes
set(PROCESS src/ProcessObservable.cpp src/DeepInelasticScattering.cpp src/F2.cpp src/FL.cpp)

# HQCDP class
set(HQCDP src/HQCDP.cpp)

add_library(hvqcd_lib ${YANGMILLS} ${HVQCD} ${SPECTRA} ${PROCESS} ${U1NNMODE} ${HQCDP} ${SCHRODINGER} ${COLNEW} ${QUADPACK})

# Fit background routines
# fitYangMills 
add_executable(fitYangMills.exe physics/fitYangMills.cpp)
target_link_libraries(fitYangMills.exe armadillo hvqcd_lib)

# Fit HVQCD Global
add_executable(fitHVQCDGlobal.exe physics/fitHVQCDGlobal.cpp)
target_link_libraries(fitHVQCDGlobal.exe armadillo hvqcd_lib)

# Fit HVQCD Global
add_executable(fitHVQCDGlobalWitha0.exe physics/fitHVQCDGlobalWitha0.cpp)
target_link_libraries(fitHVQCDGlobalWitha0.exe armadillo hvqcd_lib)

# Fit HVQCD All with f2
add_executable(fitHVQCDAllWithf2.exe physics/fitHVQCDAllWithf2.cpp)
target_link_libraries(fitHVQCDAllWithf2.exe armadillo hvqcd_lib)

# Fit HVQCD Global a1 and a2 fixed
add_executable(fitHVQCDGlobal_a1_a2_fixed.exe physics/fitHVQCDGlobal_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDGlobal_a1_a2_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD All with f2 a1 and a2 fixed
add_executable(fitHVQCDAllWithf2_a1_a2_fixed.exe physics/fitHVQCDAllWithf2_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDAllWithf2_a1_a2_fixed.exe armadillo hvqcd_lib)

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

add_executable(HVQCDRatiosWitha0.exe physics/HVQCDRatiosWitha0.cpp)
target_link_libraries(HVQCDRatiosWitha0.exe armadillo hvqcd_lib)

# Save background fields, potentials and schrodinger potentials of the fluctuations
add_executable(savePotentials.exe physics/savePotentials.cpp)
target_link_libraries(savePotentials.exe armadillo hvqcd_lib)

# Add test_schrodinger executable
add_executable(test_schrodinger.exe  tests/test_schrodinger.cpp)
target_link_libraries(test_schrodinger.exe armadillo hvqcd_lib)

# Add test_HVQCD executable
add_executable(test_HVQCD.exe  tests/test_HVQCD.cpp)
target_link_libraries(test_HVQCD.exe armadillo hvqcd_lib)

# Add test_HVQCD executable
add_executable(test_Kernel.exe  tests/test_Kernel.cpp)
target_link_libraries(test_Kernel.exe armadillo hvqcd_lib)

# Add test_exp_factor
add_executable(test_Wf.exe tests/test_Wf.cpp)
target_link_libraries(test_Wf.exe armadillo blas hvqcd_lib)

# Add test_EfficientNelderMead
add_executable(test_NelderMead.exe  tests/test_NelderMead.cpp)

# Add test_U1NNMode
add_executable(test_U1NNMode.exe tests/test_U1NNMode.cpp)
target_link_libraries(test_U1NNMode.exe hvqcd_lib blas armadillo)

# Add DIS fit
add_executable(fitDISHVQCD.exe physics/fitDISHVQCD.cpp)
target_link_libraries(fitDISHVQCD.exe armadillo blas  hvqcd_lib)

# Add HVQCD spectrum computation
add_executable(HVQCDSpectrum.exe physics/HVQCDSpectrum.cpp)
target_link_libraries(HVQCDSpectrum.exe armadillo hvqcd_lib)

# Add findReggeTrajectories computation
add_executable(findReggeTrajectories.exe physics/findReggeTrajectories.cpp)
target_link_libraries(findReggeTrajectories.exe armadillo hvqcd_lib)

# Add findGs computation
add_executable(findGs.exe physics/findGs.cpp)
target_link_libraries(findGs.exe armadillo blas hvqcd_lib)