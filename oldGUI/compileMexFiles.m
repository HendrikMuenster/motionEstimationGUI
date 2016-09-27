clear all;close all;clc;

mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L2L2OpticalFlowCPP.cpp' 'algorithms/mex/tools.cpp'
mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L2L2MassPreservationCPP.cpp' 'algorithms/mex/tools.cpp'

mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L2TVOpticalFlowCPP.cpp' 'algorithms/mex/tools.cpp'
mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L2TVMassPreservationCPP.cpp' 'algorithms/mex/tools.cpp'

mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L1TVOpticalFlowCPP.cpp' 'algorithms/mex/tools.cpp'
mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L1TVMassPreservationCPP.cpp' 'algorithms/mex/tools.cpp'

mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L1L2OpticalFlowCPP.cpp' 'algorithms/mex/tools.cpp'
mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L1L2MassPreservationCPP.cpp' 'algorithms/mex/tools.cpp'

mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp -O3 -funroll-loops" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp -O3 -funroll-loops" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L1TVL2OpticalFlowCPP.cpp' 'algorithms/mex/tools.cpp'
mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp -O3 -funroll-loops" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp -O3 -funroll-loops" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L1TVTVOpticalFlowCPP.cpp' 'algorithms/mex/tools.cpp'

mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp -O3 -funroll-loops" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp -O3 -funroll-loops" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L2TVBregOpticalFlowCPP.cpp' 'algorithms/mex/tools.cpp'
%%
mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L1TVOpticalFlowNonlinearCPP.cpp' 'algorithms/mex/tools.cpp'  'algorithms/mex/linearInterpolation.cpp'  'algorithms/mex/cubicInterpolation.cpp'