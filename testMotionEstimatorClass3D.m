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

image1 = zeros(100,100,100);
image2 = zeros(100,100,100);

image1(20:40,40:60,20:40) = 1;
image2(60:80,60:80,60:80) = 1;

cut1 = squeeze(image1(:,10,:));
cut2 = squeeze(image2(:,10,:));

figure(1);imagesc(cut1);axis image;
figure(2);imagesc(cut2);axis image;

%%
u = cat(4,image1,image2);

tol = 1e-5;
alpha = 0.01;

verbose = 2;
dataTerm = 'L1';
regularizerTerm = 'Huber';
doGradientConstancy = 0;
steplength = 0.8;
numberOfWarps = 5;

%% 

motionEstimator = motionEstimatorClass(u,tol,alpha,'dataTerm',dataTerm,'regularizerTerm',regularizerTerm,'doGradientConstancy',doGradientConstancy,'steplength',steplength,'numberOfWarps',numberOfWarps);
motionEstimator.init;

motionEstimator.verbose = verbose;

%%
tic;motionEstimator.runPyramid;toc;

%%
v = motionEstimator.getResult;

%%
for i=1:3
    for j=1:size(v,3)
        figure(1);imagesc(v(:,:,j,1,i));axis image;title(['Field ',num2str(j),', component ',num2str(i)]);colorbar;pause
    end
end
error = absoluteError(squeeze(v(:,:,1,:)),flowGT)