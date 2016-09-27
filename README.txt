% MotionEstimation GUI 
%
% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 2.0
% Date: 2016-09-27

% All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author and University of Muenster not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.

If you use this framwork, please use the following citation
@PhdThesis{Dir15,
  author       = {Dirks, H.},
  title        = {Variational Methods for Joint Motion Estimation and Image Reconstruction},
  school       = {Institute for Computational and Applied Mathematics University of Muenster},
  month        = {june},
  year         = {2015},
  keywords     = {Image Processing, Image Reconstruction, Motion Estimation, Joint Image Reconstruction and Motion Estimation, Total Variation, Optical Flow, Microscopy Imaging, Biomedical Imaging, Temporal Inpainting, Variational Methods, Primal-Dual Methods},
  url          = \{/2015/Dir15},
}

Install Instructions:
Make sure you have downloaded a copy of flexBox from www.flexbox.im and this copy is in your MATALB path. 
Moreover, you should compile the C++-module for runtime reasons! Afterwards, run example.m

Further copyright notes:
- The files "computeColor.m", "writeFlowFile.m" and the original version of "flowToColorV2.m" have been developed by Deqing Sun, Department of Computer Science, Brown University (dqsun@cs.brown.edu)
- The dataset in the "data" folder comes from the Middlebury optical flow database (http://vision.middlebury.edu/flow/)

List of changes:
## Version 2.0 ##
- Motion estimator completely rewritten:
-- Based on flexBox
-- CUDA and C++ support
-- Large-scale motion estimation based on coarse-to-fine warping
## Version 1.4 ##
- Bugfixes
## Version 1.3 ##
- Bugfixes
- Missing files added
## Version 1.2 ##
- Bugfixes
## Version 1.1 ##
- Add a coarse-to-fine warping for all algorithms and set the maxIterations to default number of 10. Coarse-to-fine warping can be turned off by setting the "Pyramid steps" 
  to 1 and increasing the "Max Iterations"
- Added a L1-TV optical flow method with nonlinear data attachment term (in fact it's still a linearization, but closer to the classical attachment term)





