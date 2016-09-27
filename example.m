%% clean up

clear all;close all;clc;

addpath(genpath(cd));

%% How to:
% make sure you have downloaded a copy of flexBox from www.flexbox.im and
% this copy is in your MATALB path. Moreover, you should compile the
% C++-module for runtime reasons!
% 
% Simply edit the following section for different motion estimation models

%% parameters
dataTerm = 'L1'; %set to L1 or L2
regularizerTerm = 'L2'; %set to L2, L1 or Huber

weight = 0.01; %regulates the weight of the regularizer: higher values generate smoother results 

tol = 1e-5; %error tolerance. Don't change

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


flowGT = readFlowFile(['data',filesep,'flow10.flo']);
flowGT(abs(flowGT)>1e2) = 0;
figure(99);clf;imagesc(flowToColorV2(cat(3,flowGT(:,:,1),flowGT(:,:,2))));

%% assemble input sequence

u = cat(3,image1,image2);

%% create algorithm and solve problem
motionEstimator = motionEstimatorClass(u,tol,weight,'dataTerm',dataTerm,'regularizerTerm',regularizerTerm);
motionEstimator.init;

motionEstimator.verbose = 2;

tic;motionEstimator.runPyramid;toc;

%% get result
v = motionEstimator.getResult; %get result from class
v = squeeze(v(:,:,1,:)); %extract flow between first and second frame and remove singleton dimension

figure(1);clf;imagesc(flowToColorV2(cat(3,v(:,:,1),v(:,:,2))));