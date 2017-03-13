% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 1.00
% Date: 2016-09-20

% MotionEstimator is copyright ©2016 by Hendrik Dirks and is distributed under the terms of the GNU General Public License (GPL) version 3 (or later).
%
% If you plan to distribute the software (commercially or not). Please contact me for more information.


clear all;close all;clc;
addpath(genpath(cd));
%%
image1 = imread(['data',filesep,'frame10.png']);
image2 = imread(['data',filesep,'frame11.png']);

if (size(image1,3)>1)
    image1 = rgb2gray(image1);
    image2 = rgb2gray(image2);
end

image1 = im2double(image1);
image2 = im2double(image2);

%resize image
image1 = imresize(image1,1);
image2 = imresize(image2,1);

figure(1);imagesc(image1)
figure(2);imagesc(image2)

flowGT = readFlowFile(['data',filesep,'flow10.flo']);
flowGT2 = flowGT;
flowGT2(abs(flowGT)>1e2) = 0;
figure(99);clf;imagesc(flowToColorV2(cat(3,flowGT2(:,:,1),flowGT2(:,:,2))));


%%
u = cat(3,image1,image2);

tol = 1e-5;
alpha = 0.01;

verbose = 2;
dataTerm = 'L1';
regularizerTerm = 'Huber';
doGradientConstancy = 0;
steplength = 0.9;
numberOfWarps = 5;

%% 

motionEstimator = motionEstimatorClass(u,tol,alpha,'verbose',verbose,'dataTerm',dataTerm,'regularizerTerm',regularizerTerm,'doGradientConstancy',doGradientConstancy,'steplength',steplength,'numberOfWarps',numberOfWarps);
motionEstimator.init;


%%
tic;motionEstimator.runPyramid;toc;

%%
v = motionEstimator.getResult;

%%
tmp = v(:,:,1,1);
v(:,:,1,1) = v(:,:,1,2);
v(:,:,1,2) = tmp;
error = absoluteError(squeeze(v(:,:,1,:)),flowGT)