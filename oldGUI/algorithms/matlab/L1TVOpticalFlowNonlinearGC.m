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
function [x,y] = L1TVOpticalFlowNonlinearGC(image1,image2,tol,lambda,varargin)
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
    
    Kx(:,5) = 0;
    Kx(:,6) = 0;
    Kx(:,7) = 0;
    
    Kty = x;

    reverseStr = [];

    for steps=1:5
        %Warp Image 1
        xReshape = reshape(x,[M,N,2]);
        
        xReshape(:,:,1) = medfilt2(xReshape(:,:,1), [5,5], 'symmetric');
        xReshape(:,:,2) = medfilt2(xReshape(:,:,2), [5,5], 'symmetric');
        
        origX = imageGrid(:,:,1);
        origY = imageGrid(:,:,2);
        
        shiftX = imageGrid(:,:,1) + xReshape(:,:,1);
        shiftY = imageGrid(:,:,2) + xReshape(:,:,2);
        
        markerOutOfGrid = (shiftX>size(shiftX,2)) + (shiftX<1) + (shiftY>size(shiftX,1)) + (shiftY<1);
        markerOutOfGrid = markerOutOfGrid > 0;
        %these are for I2(x+\tilde{v})
        idxx = max(1,min(N,shiftX));
        idxm = max(1,min(N,shiftX-0.5));
        idxp = max(1,min(N,shiftX+0.5));

        idyy = max(1,min(M,shiftY));
        idym = max(1,min(M,shiftY-0.5));
        idyp = max(1,min(M,shiftY+0.5));
        
        %these are for I1(x)
        idxx2 = max(1,min(N,origX));
        idxm2 = max(1,min(N,origX-0.5));
        idxp2 = max(1,min(N,origX+0.5));

        idyy2 = max(1,min(M,origY));
        idym2 = max(1,min(M,origY-0.5));
        idyp2 = max(1,min(M,origY+0.5));
        
        I1x = interp2(image1,idxp2,idyy2,'bicubic') - interp2(image1,idxm2,idyy2,'bicubic');
        I1y = interp2(image1,idxx2,idyp2,'bicubic') - interp2(image1,idxx2,idym2,'bicubic');
        
        I2w = interp2(image2,idxx,idyy,'bicubic');
        I2wx = interp2(image2,idxp,idyy,'bicubic') - interp2(image2,idxm,idyy,'bicubic');
        I2wy = interp2(image2,idxx,idyp,'bicubic') - interp2(image2,idxx,idym,'bicubic');
        
        %second order derivatives
        I2wxx = interp2(I2wx,idxp2,idyy2,'bicubic') - interp2(I2wx,idxm2,idyy2,'bicubic');
        I2wxy = interp2(I2wx,idxx2,idyp2,'bicubic') - interp2(I2wx,idxx2,idym2,'bicubic');
        I2wyx = interp2(I2wy,idxp2,idyy2,'bicubic') - interp2(I2wy,idxm2,idyy2,'bicubic');
        I2wyy = interp2(I2wy,idxx2,idyp2,'bicubic') - interp2(I2wy,idxx2,idym2,'bicubic');
        
        %sum up terms in 
        ut = I2w - image1 - I2wx.*xReshape(:,:,1) - I2wy.*xReshape(:,:,2);
        
        
        uxt = I2wx - I1x - I2wxx.*xReshape(:,:,1) - I2wxy.*xReshape(:,:,2);
        uyt = I2wy - I1y - I2wyx.*xReshape(:,:,1) - I2wyy.*xReshape(:,:,2);
        
        %set data term to zero
        ut(markerOutOfGrid) = 0;
        I2wx(markerOutOfGrid) = 0;
        I2wy(markerOutOfGrid) = 0;
        
        I2wxx(markerOutOfGrid) = 0;
        I2wxy(markerOutOfGrid) = 0;
        I2wyx(markerOutOfGrid) = 0;
        I2wyy(markerOutOfGrid) = 0;
        uxt(markerOutOfGrid) = 0;
        uyt(markerOutOfGrid) = 0;
        
        ut = ut(:);
        ux = I2wx(:);
        uy = I2wy(:);
        
        
        
        sigma1 = 2 / stepsize(2);
        sigma2 = max(abs(I2wxx(:))) + max(abs(I2wxy(:)));
        sigma3 = max(abs(I2wyx(:))) + max(abs(I2wyy(:)));
        sigma4 = max(abs(ux(:))) + max(abs(uy(:)));
        
        tau1 = 4 / stepsize(2) + max(abs(I2wxx(:))) + max(abs(I2wyx(:))) + max(abs(ux(:)));
        tau2 = 4 / stepsize(2) + max(abs(I2wxy(:))) + max(abs(I2wyy(:))) + max(abs(uy(:)));
        
        sigma1 = 1/sigma1;
        sigma2 = 1/sigma2;
        sigma3 = 1/sigma3;
        sigma4 = 1/sigma4;
        
        tau1 = 1/tau1;
        tau2 = 1/tau2;
        
        clear image1WarpGradient;
        
        iterations = 1;err=1;
        while err>tol && iterations < maxIt
            yOld = y;
            xOld = x;

            KtyOld = Kty;

            Kty(:,1) = Dx2D'*y(:,1) + Dy2D'*y(:,2) + I2wxx(:).*y(:,5) + I2wyx(:).*y(:,6) + ux(:).*y(:,7);
            Kty(:,2) = Dx2D'*y(:,3) + Dy2D'*y(:,4) + I2wxy(:).*y(:,5) + I2wyy(:).*y(:,6) + uy(:).*y(:,7);

            x(:,1) = x(:,1) - tau1*Kty(:,1);
            x(:,2) = x(:,2) - tau2*Kty(:,2);

            KxOld = Kx;
            Kx(:,1) = Dx2D*x(:,1);
            Kx(:,2) = Dy2D*x(:,1);
            Kx(:,3) = Dx2D*x(:,2);
            Kx(:,4) = Dy2D*x(:,2);
            
            %new data terms
            Kx(:,5) = I2wxx(:).*x(:,1) + I2wxy(:).*x(:,2);
            Kx(:,6) = I2wyx(:).*x(:,1) + I2wyy(:).*x(:,2);
            Kx(:,7) = ux.*x(:,1) + uy.*x(:,2);
            
            %dual prox
            yTilde(:,1:4)  = y(:,1:4) + sigma1*(2*Kx(:,1:4) - KxOld(:,1:4));
            yTilde(:,5)  = y(:,5) + sigma2*(2*Kx(:,5) - KxOld(:,5));
            yTilde(:,6)  = y(:,6) + sigma3*(2*Kx(:,6) - KxOld(:,6));
            yTilde(:,7)  = y(:,7) + sigma4*(2*Kx(:,7) - KxOld(:,7));
            
            par = 0;
            y(:,5) = max(-par,min(par,yTilde(:,5) + sigma2*uxt(:)));
            y(:,6) = max(-par,min(par,yTilde(:,6) + sigma3*uyt(:)));
            y(:,7) = max(-1,min(1,yTilde(:,7) + sigma4*ut(:)));
            
            %regularizers
            
            if (typeNorm == 1)
                %anisotropic
                y = max(-lambda,min(lambda,y + sigma1*(2*Kx-KxOld)));
            elseif(typeNorm == 4)
                %Huber norm
                huberFactor = 1 / (1 + sigma1* huberEpsilon / lambda);
                
                yForward = yTilde(:,1:4) * huberFactor;

                norm1 = sqrt(yForward(:,1).^2 + yForward(:,2).^2);
                norm2 = sqrt(yForward(:,3).^2 + yForward(:,4).^2);

                y(:,1:4) = yForward ./ max(1,cat(2,norm1/lambda,norm1/lambda,norm2/lambda,norm2/lambda));
            else
                %mixed isotropic
                yForward = y + sigma*(2*Kx-KxOld);

                norm1 = sqrt(yForward(:,1).^2 + yForward(:,2).^2);
                norm2 = sqrt(yForward(:,3).^2 + yForward(:,4).^2);

                y = yForward ./ max(1,cat(2,norm1/lambda,norm1/lambda,norm2/lambda,norm2/lambda));
            end
            
            if (mod(iterations,100)==1)
                xDiff = xOld - x;
                yDiff = yOld - y;

                %primal residual
                p1 = xDiff(:,1)/tau1 - (KtyOld(:,1)-Kty(:,1));
                p2 = xDiff(:,2)/tau2 - (KtyOld(:,2)-Kty(:,2));
                %dual residual
                d1 = yDiff(:,1:4)/sigma1 - (KxOld(:,1:4)-Kx(:,1:4));
                d2 = yDiff(:,5)/sigma2 - (KxOld(:,5)-Kx(:,5));
                d3 = yDiff(:,6)/sigma3 - (KxOld(:,6)-Kx(:,6));
                d4 = yDiff(:,7)/sigma4 - (KxOld(:,7)-Kx(:,7));


                p=sum(abs(p1(:))) + sum(abs(p2(:)));
                d=sum(abs(d1(:))) + sum(abs(d2(:))) + sum(abs(d3(:))) + sum(abs(d4(:)));

                err = (p+d) / nPx;
                
                reverseStr = printToCmd( reverseStr,sprintf('L1-TV Large-Scale Optical Flow, Iteration: #%d : Residual %f',iterations,err) );
                
                %pause
            end

            iterations = iterations + 1;
        end
    end
    
    %clear comand output
    printToCmd( reverseStr,'');
    
    x = reshape(x,[size(image1),2]);
    y = reshape(y,[size(image1),7]);
        
end

