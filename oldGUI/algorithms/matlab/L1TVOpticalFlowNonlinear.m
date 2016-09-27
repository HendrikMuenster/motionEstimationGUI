% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 1.0
% Date: 2015-06-17

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
%
function [x,y] = L1TVOpticalFlowNonlinear(image1,image2,tol,lambda,varargin)
    %insert is a list of images, output a list of velocity fields
    %returns primal and dual variable
    vararginParser;
	nPx = numel(image1);
    
    %allocate grid for warping
    M = size(image1,1);
    N = size(image1,2);
    [imageGrid(:,:,1),imageGrid(:,:,2)] = meshgrid(1:N,1:M);

	D2Dp = generateForwardGradientND( size(image1),[1,1]);
	
	Dx2D = D2Dp(0*nPx+1:1*nPx,:);
	Dy2D = D2Dp(1*nPx+1:2*nPx,:);
    
    Dx2D = 1/stepsize(2) * Dx2D;
    Dy2D = 1/stepsize(2) * Dy2D;
    
    nPx = numel(image1);
    tau = stepsize(2)/4;
    sigma = stepsize(2)/2;

    if (~exist('x','var'))
        x = zeros(nPx,2);
    end
    if (~exist('y','var'))
        y = zeros(nPx,4);
    end
    if (~exist('maxIt','var'))
        maxIt = 50000;
    end
    if (~exist('typeNorm','var'))
        typeNorm = 2;
    end
    
    Kx(:,1) = Dx2D*x(:,1);
    Kx(:,2) = Dy2D*x(:,1);
    Kx(:,3) = Dx2D*x(:,2);
    Kx(:,4) = Dy2D*x(:,2);
    
    Kty = x;

    reverseStr = [];

    for steps=1:5
        %Warp Image 1
        xReshape = reshape(x,[M,N,2]);
        
        origX = imageGrid(:,:,1);
        origY = imageGrid(:,:,2);
        
        shiftX = imageGrid(:,:,1) + xReshape(:,:,1);
        shiftY = imageGrid(:,:,2) + xReshape(:,:,2);

        idxx = max(1,min(N,shiftX));
        idxm = max(1,min(N,shiftX-0.5));
        idxp = max(1,min(N,shiftX+0.5));

        idyy = max(1,min(M,shiftY));
        idym = max(1,min(M,shiftY-0.5));
        idyp = max(1,min(M,shiftY+0.5));
        
        idxx2 = max(1,min(N,origX));
        idxm2 = max(1,min(N,origX-0.5));
        idxp2 = max(1,min(N,origX+0.5));

        idyy2 = max(1,min(M,origY));
        idym2 = max(1,min(M,origY-0.5));
        idyp2 = max(1,min(M,origY+0.5));

        image2Warp = interp2(image2,idxx,idyy,'bicubic');
        image1WarpGradient(:,:,1) = interp2(image2,idxp,idyy,'bicubic') - interp2(image2,idxm,idyy,'bicubic');
        image1WarpGradient(:,:,2) = interp2(image2,idxx,idyp,'bicubic') - interp2(image2,idxx,idym,'bicubic');
        image1WarpGradient(:,:,3) = image2Warp - image1;
        
        image2Gradient(:,:,1) = interp2(image1,idxp2,idyy2,'bicubic') - interp2(image1,idxm2,idyy2,'bicubic');
        image2Gradient(:,:,2) = interp2(image1,idxx2,idyp2,'bicubic') - interp2(image1,idxx2,idym2,'bicubic');
        

        image1WarpGradient(1,:,:) = 0;
        image1WarpGradient(end,:,:) = 0;
        image1WarpGradient(:,1,:) = 0;
        image1WarpGradient(:,end,:) = 0;
        
        image2Gradient(1,:,:) = 0;
        image2Gradient(end,:,:) = 0;
        image2Gradient(:,1,:) = 0;
        image2Gradient(:,end,:) = 0;
        
        %image1WarpGradient(:,:,1) = 0.5 * (image1WarpGradient(:,:,1) + image2Gradient(:,:,1));
        %image1WarpGradient(:,:,2) = 0.5 * (image1WarpGradient(:,:,2) + image2Gradient(:,:,2));
        
        %mask = fspecial('gaussian', [7, 7], 0.5);
        
        %image1WarpGradient(:,:,1) = imfilter(image1WarpGradient(:,:,1),mask,'replicate');
        %image1WarpGradient(:,:,2) = imfilter(image1WarpGradient(:,:,2),mask,'replicate');
        
        
        ut = image1WarpGradient(:,:,3) - image1WarpGradient(:,:,1).*xReshape(:,:,1) - image1WarpGradient(:,:,2).*xReshape(:,:,2);
        ut = ut(:);
        
        ux = image1WarpGradient(:,:,1);ux = ux(:);
        uy = image1WarpGradient(:,:,2);uy = uy(:);
        
        betaQuad = max(1e-9,ux.^2+uy.^2);
        tauBetaQuad = tau*betaQuad;

        teiler(:,1) = ux ./ betaQuad;
        teiler(:,2) = uy ./ betaQuad;

        tauGrad(:,1) = tau*ux;
        tauGrad(:,2) = tau*uy;
        
        clear image1WarpGradient;
        
        iterations = 1;err=1;
        while err>tol && iterations < maxIt
            yOld = y;
            xOld = x;

            if (mod(iterations,100)==1)
                reverseStr = printToCmd( reverseStr,sprintf('L1-TV Large-Scale Optical Flow, Iteration: #%d : Residual %f',iterations,err) );
                %pause
            end

            KtyOld = Kty;

            Kty(:,1) = Dx2D'*y(:,1) + Dy2D'*y(:,2);
            Kty(:,2) = Dx2D'*y(:,3) + Dy2D'*y(:,4);

            xHat = xOld - tau*Kty;

            rho = ut + ux.*xHat(:,1) + uy.*xHat(:,2);

            c1 = rho<= -tauBetaQuad;
            c2 = rho>= tauBetaQuad;

            c1c2 = c1-c2;
            c3 = 1-c1-c2;

            %primal prox
            x = xHat + tauGrad .* repmat(c1c2,1,2) -  teiler .* repmat(rho.*c3,1,2);

            KxOld = Kx;
            Kx(:,1) = Dx2D*x(:,1);
            Kx(:,2) = Dy2D*x(:,1);
            Kx(:,3) = Dx2D*x(:,2);
            Kx(:,4) = Dy2D*x(:,2);

            %dual prox
            if (typeNorm == 1)
                %anisotropic
                y = max(-lambda,min(lambda,y + sigma*(2*Kx-KxOld)));
            elseif(typeNorm == 4)
                %Huber norm
                huberFactor = 1 / (1 + sigma* huberEpsilon / lambda);
                
                yForward = (y + sigma*(2*Kx-KxOld)) * huberFactor;

                norm1 = sqrt(yForward(:,1).^2 + yForward(:,2).^2);
                norm2 = sqrt(yForward(:,3).^2 + yForward(:,4).^2);

                y = yForward ./ max(1,cat(2,norm1/lambda,norm1/lambda,norm2/lambda,norm2/lambda));
            else
                %mixed isotropic
                yForward = y + sigma*(2*Kx-KxOld);

                norm1 = sqrt(yForward(:,1).^2 + yForward(:,2).^2);
                norm2 = sqrt(yForward(:,3).^2 + yForward(:,4).^2);

                y = yForward ./ max(1,cat(2,norm1/lambda,norm1/lambda,norm2/lambda,norm2/lambda));
                
%                 %fully isotropic
%                 yForward = y + sigma*(2*Kx-KxOld);
% 
%                 norm = sqrt(yForward(:,1).^2 + yForward(:,2).^2 + yForward(:,3).^2 + yForward(:,4).^2 + 1e-8);
% 
%                 y = yForward ./ max(1,repmat(norm/lambda,1,4));
            end

            xDiff = xOld - x;
            yDiff = yOld - y;

            %primal residual
            p = xDiff/tau - (KtyOld-Kty);
            %dual residual
            d = yDiff/sigma - (KxOld-Kx);

            p=sum(abs(p(:)));
            d=sum(abs(d(:)));

            err = (p+d) / nPx;

            iterations = iterations + 1;
        end
    end
    
    %clear comand output
    printToCmd( reverseStr,'');
    
    x = reshape(x,[size(image1),2]);
    y = reshape(y,[size(image1),4]);
        
end

