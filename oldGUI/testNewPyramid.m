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


flowGT = readFlowFile(['data',filesep,'flow10.flo']);
flowGT2 = flowGT;
flowGT2(abs(flowGT)>1e2) = 0;
figure(99);clf;imagesc(flowToColorV2(cat(3,flowGT2(:,:,1),flowGT2(:,:,2))));



u = cat(3,image1,image2);

tol = 1e-5;

%u = u(100:300,100:300,:);
%flowGT = flowGT(100:300,100:300,:);

dimsU = size(u);

%%

doGaussianSmoothing = 1;
medianFiltering = 1;
medianFilteringTerm = 0;
adjustStepsize = 0;
gradientConstancy = 0;

alpha = 0.015;
numDualVars = 5;

%%
tic;
%[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1TVOpticalFlowNonlinearCuda',numDualVars,'useCPP',1,'maxIt',5000,'numsteps',1000,'adjustStepsize',adjustStepsize,'medianFiltering',medianFiltering,'doGaussianSmoothing',doGaussianSmoothing,'steplength',0.8,'typeNorm',4,'huberEpsilon',0.01);
%[x,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1TVOpticalFlowNonlinear',numDualVars,'useCPP',1,'gradientConstancy',gradientConstancy,'numsteps',1000,'adjustStepsize',adjustStepsize,'medianFiltering',medianFiltering,'doGaussianSmoothing',doGaussianSmoothing,'steplength',0.8,'typeNorm',4,'huberEpsilon',0.01);
%[x,~] = motionEstimationPyramid(u,dimsU,1e-4,alpha,'L1TVOpticalFlowNonlinear',numDualVars,'gradientConstancy',1);
[x,~] = motionEstimationPyramidFlexBox(u,dimsU,tol,alpha,'L1TVOpticalFlowNonlinear',numDualVars,'useCPP',1,'doGradientConstancy',gradientConstancy,'adjustStepsize',adjustStepsize,'medianFiltering',medianFiltering,'medianFilteringTerm',medianFilteringTerm,'doGaussianSmoothing',doGaussianSmoothing);
toc;

%%


motionEstimator = motionEstimatorClass(u,tol,alpha);
motionEstimator.init;
motionEstimator.resetImages(u);

%%
motionEstimator.runPyramid
%%
x = motionEstimator.getResult;

%%
alpha2 = 0.05;
tic;
%profile on;
[x2,~] = motionEstimationPyramid(u,dimsU,tol,alpha2,'L1TVOpticalFlowNonlinear',7,'adjustStepsize',adjustStepsize,'medianFiltering',medianFiltering,'doGaussianSmoothing',doGaussianSmoothing);
%profile off;
%profile viewer;
toc;

% image1s = imresize(u(:,:,1),0.2);
% image2s = imresize(u(:,:,2),0.2);
% u = cat(3,image1s,image2s);
% dimsU = size(u);
% [x2,~] = motionEstimationPyramid(u,dimsU,tol,alpha,'L1TVOpticalFlowNonlinear',numDualVars,'useCPP',1,'maxIt',200,'numsteps',100,'adjustStepsize',adjustStepsize,'medianFiltering',medianFiltering,'doGaussianSmoothing',doGaussianSmoothing,'steplength',0.95,'typeNorm',4,'huberEpsilon',0.01);


%%
% clear field
% clear field2
% 
% field = squeeze(x(:,:,1,:));
% 
% sF = size(field);
% 
% 
% field2(:,:,1) = imresize(x2(:,:,1,1),sF(1:2));
% field2(:,:,2) = imresize(x2(:,:,1,2),sF(1:2));
% 
% err = abs(field-field2);
% err = sum(err(:)) / numel(err(:))
%%
field = squeeze(x(:,:,1,:));

figure(9);clf;imagesc(flowToColorV2(field));axis image;
figure(10);clf;imagesc(flowToColorV2(flowGT));axis image;

bord = 1;
absoluteError(field(bord:end-bord,bord:end-bord,:),flowGT(bord:end-bord,bord:end-bord,:))
angularError(field(bord:end-bord,bord:end-bord,:),flowGT(bord:end-bord,bord:end-bord,:))


%%
field = squeeze(x(:,:,1,:));

flowGT2 = flowGT;
            UNKNOWN_FLOW_THRESH = 1e9; 
            flowGT2 (flowGT2>UNKNOWN_FLOW_THRESH) = NaN;

[aae stdae aepe] = flowAngErr(flowGT2(:,:,1), flowGT2(:,:,2), field(:,:,1), field(:,:,2), 0)



% %% do some parameter tests
% 
% doGaussianSmoothing = 1;
% medianFiltering = 1;
% adjustStepsize = 1;
% 
% num = 1;
% 
% for alpha=[0.005,0.01,0.015,0.02]
%     [x,~] = motionEstimationPyramidFlexBox(u,dimsU,tol,alpha,'L1TVOpticalFlowNonlinear',4,'useCPP',1,'maxIt',2000,'numsteps',1000,'adjustStepsize',adjustStepsize,'medianFiltering',medianFiltering,'doGaussianSmoothing',doGaussianSmoothing,'steplength',0.8,'typeNorm',4,'huberEpsilon',0.01);
%     [aae stdae aepe] = flowAngErr(flowGT2(:,:,1), flowGT2(:,:,2), field(:,:,1), field(:,:,2), 0);
%     
%     listaepe(num) = aepe;
%     
%     num = num + 1;
% end