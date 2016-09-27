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

%image1 = image1(1:200,1:200);
%image2 = image2(1:200,1:200);

flowGT = readFlowFile(['data',filesep,'flow10.flo']);
flowGT2 = flowGT;
flowGT2(abs(flowGT)>1e2) = 0;
figure(99);clf;imagesc(flowToColorV2(cat(3,flowGT2(:,:,1),flowGT2(:,:,2))));


%%
u = cat(3,image1,image2);

tol = 1e-5;
alpha = 0.05;

%% 

motionEstimator = motionEstimatorClass(u,tol,alpha);
motionEstimator.init;

motionEstimator.verbose = 1;

%%
motionEstimator.runPyramid

%%

motionEstimator.resetImages(u);

%%
motionEstimator.runLevel(1);
%%
v = motionEstimator.getResult;