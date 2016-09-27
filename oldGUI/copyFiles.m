clear all;close all;clc;

%Algorithms
copyfile('../../Dissertation/code/variationalModels/L2L2OpticalFlow.m','algorithms/matlab/L2L2OpticalFlow.m');
copyfile('../../Dissertation/code/variationalModels/L2L2MassPreservation.m','algorithms/matlab/L2L2MassPreservation.m');
copyfile('../../Dissertation/code/variationalModels/L2TVOpticalFlow.m','algorithms/matlab/L2TVOpticalFlow.m');
copyfile('../../Dissertation/code/variationalModels/L2TVMassPreservation.m','algorithms/matlab/L2TVMassPreservation.m');
copyfile('../../Dissertation/code/variationalModels/L1TVOpticalFlow.m','algorithms/matlab/L1TVOpticalFlow.m');
copyfile('../../Dissertation/code/variationalModels/L1TVMassPreservation.m','algorithms/matlab/L1TVMassPreservation.m');
copyfile('../../Dissertation/code/variationalModels/L1L2OpticalFlow.m','algorithms/matlab/L1L2OpticalFlow.m');
copyfile('../../Dissertation/code/variationalModels/L1L2MassPreservation.m','algorithms/matlab/L1L2MassPreservation.m');
copyfile('../../Dissertation/code/variationalModels/L1TVOpticalFlowNonlinear.m','algorithms/matlab/L1TVOpticalFlowNonlinear.m');

copyfile('../../Dissertation/code/variationalModels/L1TVL2OpticalFlow.m','algorithms/matlab/L1TVL2OpticalFlow.m');
copyfile('../../Dissertation/code/variationalModels/L1TVTVOpticalFlow.m','algorithms/matlab/L1TVTVOpticalFlow.m');

%Operators
copyfile('../../Dissertation/code/variationalModels/operators/generateForwardGradient2D.m','operators/generateForwardGradient2D.m');
copyfile('../../Dissertation/code/variationalModels/operators/generateCentralGradient2D.m','operators/generateCentralGradient2D.m');
copyfile('../../Dissertation/code/variationalModels/operators/generateForwardGradientND.m','operators/generateForwardGradientND.m');
copyfile('../../Dissertation/code/variationalModels/operators/generateCentralGradientND.m','operators/generateCentralGradientND.m');

%Mex Files
copyfile('../../Dissertation/code/mex/mexTools/tools.cpp','algorithms/mex/tools.cpp');
copyfile('../../Dissertation/code/mex/mexTools/tools.h','algorithms/mex/tools.h');

copyfile('../../Dissertation/code/mex/L2L2OpticalFlowSource/functions.cpp','algorithms/mex/L2L2OpticalFlowCPP.cpp');
copyfile('../../Dissertation/code/mex/L2L2MassPreservationSource/functions.cpp','algorithms/mex/L2L2MassPreservationCPP.cpp');
copyfile('../../Dissertation/code/mex/L2TVOpticalFlowSource/functions.cpp','algorithms/mex/L2TVOpticalFlowCPP.cpp');
copyfile('../../Dissertation/code/mex/L2TVMassPreservationSource/functions.cpp','algorithms/mex/L2TVMassPreservationCPP.cpp');
copyfile('../../Dissertation/code/mex/L1TVOpticalFlowSource/functions.cpp','algorithms/mex/L1TVOpticalFlowCPP.cpp');
copyfile('../../Dissertation/code/mex/L1TVMassPreservationSource/functions.cpp','algorithms/mex/L1TVMassPreservationCPP.cpp');
copyfile('../../Dissertation/code/mex/L1L2OpticalFlowSource/functions.cpp','algorithms/mex/L1L2OpticalFlowCPP.cpp');
copyfile('../../Dissertation/code/mex/L1L2MassPreservationSource/functions.cpp','algorithms/mex/L1L2MassPreservationCPP.cpp');
%Breg
copyfile('../../Dissertation/code/mex/L2TVBregOpticalFlowSource/functions.cpp','algorithms/mex/L2TVBregOpticalFlowCPP.cpp');
%2nd Reg
copyfile('../../Dissertation/code/mex/L1TVL2OpticalFlowSource/functions.cpp','algorithms/mex/L1TVL2OpticalFlowCPP.cpp');
copyfile('../../Dissertation/code/mex/L1TVTVOpticalFlowSource/functions.cpp','algorithms/mex/L1TVTVOpticalFlowCPP.cpp');

copyfile('../../Dissertation/code/mex/L1TVOpticalFlowNonlinearSource/functions.cpp','algorithms/mex/L1TVOpticalFlowNonlinearCPP.cpp');
copyfile('../../Dissertation/code/mex/mexTools/linearInterpolation.h','algorithms/mex/linearInterpolation.h');
copyfile('../../Dissertation/code/mex/mexTools/linearInterpolation.cpp','algorithms/mex/linearInterpolation.cpp');
copyfile('../../Dissertation/code/mex/mexTools/cubicInterpolation.h','algorithms/mex/cubicInterpolation.h');
copyfile('../../Dissertation/code/mex/mexTools/cubicInterpolation.cpp','algorithms/mex/cubicInterpolation.cpp');


%Core pyramid file
copyfile('../../Dissertation/code/variationalModels/motionEstimationPyramid.m','motionEstimationPyramid.m');

%Tools
copyfile('../../matlabToolsPath/flowToColorV2.m','tools/flowToColorV2.m');
copyfile('../../matlabToolsPath/computeColor.m','tools/computeColor.m');
copyfile('../../matlabToolsPath/writeFlowFile.m','tools/writeFlowFile.m');
copyfile('../../matlabToolsPath/readFlowFile.m','tools/readFlowFile.m');
copyfile('../../Dissertation/code/variationalModels/vararginParser.m','tools/vararginParser.m');
copyfile('../../matlabToolsPath/initVar.m','tools/initVar.m');