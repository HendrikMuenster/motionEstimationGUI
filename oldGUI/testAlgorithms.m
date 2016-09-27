clear all;close all;clc;

addpath(genpath(cd));

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

u = cat(3,image1,image2);
dimsU = size(u);
tol = 1e-5;


%%
alpha = 0.1;
numDualVars = 4;
tic;
[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L2L2OpticalFlow',numDualVars,'maxIt',10000,'numsteps',1,'useCPP',1);
toc
figure(1);clf;imagesc(flowToColorV2(cat(3,x(:,:,1),x(:,:,2))));
%%
alpha = 0.1;
numDualVars = 5;
[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L2L2MassPreservation',numDualVars,'maxIt',10000,'numsteps',100);

figure(2);clf;imagesc(flowToColorV2(cat(3,x(:,:,1),x(:,:,2))));
%%
alpha = 0.0005;
numDualVars = 4;
[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L2TVOpticalFlow',numDualVars,'maxIt',10000,'numsteps',100);

figure(3);clf;imagesc(flowToColorV2(cat(3,x(:,:,1),x(:,:,2))));
%%
alpha = 0.0003;
numDualVars = 5;
[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L2TVMassPreservation',numDualVars,'maxIt',1000,'numsteps',100);

figure(4);clf;imagesc(flowToColorV2(cat(3,x(:,:,1),x(:,:,2))));
%%
alpha = 0.01;
numDualVars = 4;
[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1TVOpticalFlow',numDualVars,'maxIt',1000,'numsteps',100);

figure(5);clf;imagesc(flowToColorV2(cat(3,x(:,:,1),x(:,:,2))));
%%
alpha = 0.01;
numDualVars = 5;
[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1TVMassPreservation',numDualVars,'maxIt',1000,'numsteps',100);

figure(6);clf;imagesc(flowToColorV2(cat(3,x(:,:,1),x(:,:,2))));
%%
alpha = 0.2;
numDualVars = 4;
[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1L2OpticalFlow',numDualVars,'maxIt',1000,'numsteps',100);

figure(7);clf;imagesc(flowToColorV2(cat(3,x(:,:,1),x(:,:,2))));
%%
alpha = 0.2;
numDualVars = 5;
[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1L2MassPreservation',numDualVars,'maxIt',1000,'numsteps',100);

figure(8);clf;imagesc(flowToColorV2(cat(3,x(:,:,1),x(:,:,2))));
%%

u = cat(3,image1,image2);
dimsU = size(u);
tol = 1e-5;

alpha = 0.05;
numDualVars = 7;
[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1TVOpticalFlowNonlinear',numDualVars,'maxIt',1000,'numsteps',1);

figure(9);clf;imagesc(flowToColorV2(cat(3,x(:,:,1,1),x(:,:,1,2))));


