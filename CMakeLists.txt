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
set(SPECTRA src/Spectra.cpp src/Reggeon.cpp src/Kernel.cpp src/GluonKernel.cpp src/MesonKernel.cpp)

# Process and its derived classes
set(PROCESS src/Process.cpp src/Sigma.cpp src/SigmaProtonProton.cpp src/DeepInelasticScattering.cpp src/F2.cpp src/F2Photon.cpp src/FL.cpp src/SigmaGammaProton.cpp src/SigmaGammaGamma.cpp)

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

# Fit HVQCD Global With tachyon normalized
add_executable(fitHVQCDGlobalTachyonNormalized.exe physics/fitHVQCDGlobalTachyonNormalized.cpp)
target_link_libraries(fitHVQCDGlobalTachyonNormalized.exe armadillo hvqcd_lib)

# Fit HVQCD Global With a0
add_executable(fitHVQCDGlobalWitha0.exe physics/fitHVQCDGlobalWitha0.cpp)
target_link_libraries(fitHVQCDGlobalWitha0.exe armadillo hvqcd_lib)

# Fit HVQCD Global With a0 and with a1 and a2 fixed
add_executable(fitHVQCDGlobalWitha0_a1_a2_fixed.exe physics/fitHVQCDGlobalWitha0_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDGlobalWitha0_a1_a2_fixed.exe armadillo hvqcd_lib)


# Fit HVQCD Global With a0 Tachyon Normalized
add_executable(fitHVQCDGlobalWitha0TachyonNormalized.exe physics/fitHVQCDGlobalWitha0TachyonNormalized.cpp)
target_link_libraries(fitHVQCDGlobalWitha0TachyonNormalized.exe armadillo hvqcd_lib)


# Fit HVQCD Global a1 and a2 fixed
add_executable(fitHVQCDGlobal_a1_a2_fixed.exe physics/fitHVQCDGlobal_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDGlobal_a1_a2_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Global a1 = 0 and a2 = 1 fixed
add_executable(fitHVQCDGlobal_a1_0_a2_1_fixed.exe physics/fitHVQCDGlobal_a1_0_a2_1_fixed.cpp)
target_link_libraries(fitHVQCDGlobal_a1_0_a2_1_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Global with a0, a1 = 0 and a2 = 1 fixed
add_executable(fitHVQCDGlobalWitha0_a1_0_a2_1_fixed.exe physics/fitHVQCDGlobalWitha0_a1_0_a2_1_fixed.cpp)
target_link_libraries(fitHVQCDGlobalWitha0_a1_0_a2_1_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Global sc and VgIR fixed
add_executable(fitHVQCDGlobal_sc_VgIR_fixed.exe physics/fitHVQCDGlobal_sc_VgIR_fixed.cpp)
target_link_libraries(fitHVQCDGlobal_sc_VgIR_fixed.exe armadillo hvqcd_lib)

# Fit HVQCD Global sc, VgIR, a1 and a2 fixed
add_executable(fitHVQCDGlobal_sc_VgIR_a1_a2_fixed.exe physics/fitHVQCDGlobal_sc_VgIR_a1_a2_fixed.cpp)
target_link_libraries(fitHVQCDGlobal_sc_VgIR_a1_a2_fixed.exe armadillo hvqcd_lib)

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

# Add DIS fit
add_executable(fitDISHVQCDKernelFixed.exe physics/fitDISHVQCDKernelFixed.cpp)
target_link_libraries(fitDISHVQCDKernelFixed.exe armadillo blas  hvqcd_lib)

# Add DIS fit using onlye the gluon kernel
add_executable(fitDISHVQCDGluonKernelOnly.exe physics/fitDISHVQCDGluonKernelOnly.cpp)
target_link_libraries(fitDISHVQCDGluonKernelOnly.exe armadillo blas  hvqcd_lib)

# Add HVQCD spectrum computation
add_executable(HVQCDSpectrum.exe physics/HVQCDSpectrum.cpp)
target_link_libraries(HVQCDSpectrum.exe armadillo hvqcd_lib)

# Add findReggeTrajectories computation
add_executable(findReggeTrajectories.exe physics/findReggeTrajectories.cpp)
target_link_libraries(findReggeTrajectories.exe armadillo hvqcd_lib)

# Add fitSigmaProtonProton
add_executable(fitSigmaProtonProton.exe physics/fitSigmaProtonProton.cpp)
target_link_libraries(fitSigmaProtonProton.exe armadillo hvqcd_lib)

# Add fitSigmaProtonProtonKernelFixed
add_executable(fitSigmaProtonProtonKernelFixed.exe physics/fitSigmaProtonProtonKernelFixed.cpp)
target_link_libraries(fitSigmaProtonProtonKernelFixed.exe armadillo hvqcd_lib)

# Add fitSigmaGammaProtonKernelFixed
add_executable(fitSigmaGammaProtonKernelFixed.exe physics/fitSigmaGammaProtonKernelFixed.cpp)
target_link_libraries(fitSigmaGammaProtonKernelFixed.exe armadillo hvqcd_lib)

# Add fitSigmaGammaGammaKernelFixed
add_executable(fitSigmaGammaGammaKernelFixed.exe physics/fitSigmaGammaGammaKernelFixed.cpp)
target_link_libraries(fitSigmaGammaGammaKernelFixed.exe armadillo hvqcd_lib)

# Add fitSigmaGammaProton
add_executable(fitSigmaGammaProton.exe physics/fitSigmaGammaProton.cpp)
target_link_libraries(fitSigmaGammaProton.exe armadillo hvqcd_lib)

# Add fitSigmaGammaGamma
add_executable(fitSigmaGammaGamma.exe physics/fitSigmaGammaGamma.cpp)
target_link_libraries(fitSigmaGammaGamma.exe armadillo hvqcd_lib)

# Add fitAllSigma
add_executable(fitAllSigma.exe physics/fitAllSigma.cpp)
target_link_libraries(fitAllSigma.exe armadillo hvqcd_lib)

# Add fitAllSigmaKernelFixed
add_executable(fitAllSigmaFixedKernel.exe physics/fitAllSigmaFixedKernel.cpp)
target_link_libraries(fitAllSigmaFixedKernel.exe armadillo hvqcd_lib)

# Add findGs computation
add_executable(findGs.exe physics/findGs.cpp)
target_link_libraries(findGs.exe armadillo blas hvqcd_lib)

# Add findKernelIntercepts computation
add_executable(findKernelIntercepts.exe physics/findKernelIntercepts.cpp)
target_link_libraries(findKernelIntercepts.exe armadillo blas hvqcd_lib)

# test_BckPotFactor
add_executable(test_BckPotFactor.exe tests/test_BckPotFactor.cpp)
target_link_libraries(test_BckPotFactor.exe armadillo blas hvqcd_lib)

# test_SigmaProtonProton
add_executable(test_SigmaProtonProton.exe tests/test_SigmaProtonProton.cpp)
target_link_libraries(test_SigmaProtonProton.exe armadillo blas hvqcd_lib)

add_executable(predictSigmaGammaGamma.exe physics/predictSigmaGammaGamma.cpp)
target_link_libraries(predictSigmaGammaGamma.exe armadillo blas hvqcd_lib)

add_executable(predictSigmaGammaProton.exe physics/predictSigmaGammaProton.cpp)
target_link_libraries(predictSigmaGammaProton.exe armadillo blas hvqcd_lib)

add_executable(predictSigmaProtonProton.exe physics/predictSigmaProtonProton.cpp)
target_link_libraries(predictSigmaProtonProton.exe armadillo blas hvqcd_lib)

add_executable(test_gauss_kronrod_boost.exe tests/test_gauss_kronrod_boost.cpp)

add_executable(test_AdS_colnew.exe tests/test_AdS_colnew.cpp)
target_link_libraries(test_AdS_colnew.exe blas hvqcd_lib)

add_executable(fitFLSigmaGammaProtonKernelFixed.exe physics/fitFLSigmaGammaProtonKernelFixed.cpp)
target_link_libraries(fitFLSigmaGammaProtonKernelFixed.exe armadillo blas hvqcd_lib)

add_executable(fitFLSigmaGammaProton.exe physics/fitFLSigmaGammaProton.cpp)
target_link_libraries(fitFLSigmaGammaProton.exe armadillo blas hvqcd_lib)

add_executable(fitFLSigmaGammaProtonGluonKernelFixed.exe physics/fitFLSigmaGammaProtonGluonKernelFixed.cpp)
target_link_libraries(fitFLSigmaGammaProtonGluonKernelFixed.exe armadillo blas hvqcd_lib)

add_executable(fitFLSigmaGammaProtonGluonKernel.exe physics/fitFLSigmaGammaProtonGluonKernel.cpp)
target_link_libraries(fitFLSigmaGammaProtonGluonKernel.exe armadillo blas hvqcd_lib)

add_executable(fitGammaProtonGluonKernel.exe physics/fitGammaProtonGluonKernel.cpp)
target_link_libraries(fitGammaProtonGluonKernel.exe armadillo blas hvqcd_lib)

add_executable(fitGammaProtonGluonKernelFixed.exe physics/fitGammaProtonGluonKernelFixed.cpp)
target_link_libraries(fitGammaProtonGluonKernelFixed.exe armadillo blas hvqcd_lib)

add_executable(fitF2GluonKernel.exe physics/fitF2GluonKernel.cpp)
target_link_libraries(fitF2GluonKernel.exe armadillo blas hvqcd_lib)

add_executable(fitF2GluonKernelFixed.exe physics/fitF2GluonKernelFixed.cpp)
target_link_libraries(fitF2GluonKernelFixed.exe armadillo blas hvqcd_lib)

add_executable(predictStructureFunctions.exe physics/predictStructureFunctions.cpp)
target_link_libraries(predictStructureFunctions.exe armadillo blas hvqcd_lib)

add_executable(fitF2.exe physics/fitF2.cpp)
target_link_libraries(fitF2.exe armadillo blas hvqcd_lib)