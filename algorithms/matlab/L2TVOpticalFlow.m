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
function [x,y] = L2TVOpticalFlow(image1,image2,tol,lambda,varargin)
    %insert is a list of images, output a list of velocity fields
    %returns primal and dual variable
    vararginParser;
    
    [ Dx2Dc,Dy2Dc ] = generateCentralGradient2D( size(image1),[1 1]);
    
    uy = Dx2Dc*image1(:);
    ux = Dy2Dc*image1(:);
    ut = image2(:)-image1(:);
    
    clear Dx2Dc Dy2Dc;
    
    [ Dx2D,Dy2D ] = generateForwardGradient2D( size(image1),[1 1]);
    
    nPx = numel(image1);
    tau = 0.25;
    sigma = 0.5;
    
    uxut = ux.*ut;
    uyut = uy.*ut;
    
    c1 = 1+tau*ux.^2;
    c2 = tau*ux.*uy;
    c3 = 1+tau*uy.^2;
    
    teiler = c1.*c3-c2.^2+1e-8;
    
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
    
    Kx = y;
    Kty = x;

    reverseStr = [];
    
    iterations = 1;err=1;
    while err>tol && iterations < maxIt
        yOld = y;
        xOld = x;
        
        if (mod(iterations,100)==1)
            reverseStr = printToCmd( reverseStr,sprintf('L2-TV Optical Flow, Iteration: #%d : Residual %f',iterations,err) );
            
            %xR = reshape(x,[size(image1),2]);
            %figure(4);clf;imagesc(colourfulOrientationPlot(xR(:,:,1),xR(:,:,2),5));drawnow
        end
        
        KtyOld = Kty;
        
        Kty(:,1) = Dx2D'*y(:,1) + Dy2D'*y(:,2);
        Kty(:,2) = Dx2D'*y(:,3) + Dy2D'*y(:,4);
        
        xHat = x - tau*Kty;

        b1 = xHat(:,1) - tau*uxut;
        b2 = xHat(:,2) - tau*uyut;
        
        x(:,1) =  (b1.*c3-c2.*b2)./teiler;
        x(:,2) =  (b2.*c1-c2.*b1)./teiler;
        
        KxOld = Kx;
        Kx(:,1) = Dx2D*x(:,1);
        Kx(:,2) = Dy2D*x(:,1);
        Kx(:,3) = Dx2D*x(:,2);
        Kx(:,4) = Dy2D*x(:,2);
        
        %dual prox
        if (typeNorm == 1)
            %anisotropic
            y = max(-lambda,min(lambda,y + sigma*(2*Kx-KxOld)));
        elseif(typeNorm == 2)
            %mixed isotropic
            yForward = y + sigma*(2*Kx-KxOld);
            
            norm1 = sqrt(yForward(:,1).^2 + yForward(:,2).^2);
            norm2 = sqrt(yForward(:,3).^2 + yForward(:,4).^2);
            
            y = yForward ./ max(1,cat(2,norm1/lambda,norm1/lambda,norm2/lambda,norm2/lambda));
        else
            %fully isotropic
            yForward = y + sigma*(2*Kx-KxOld);
            
            norm = sqrt(yForward(:,1).^2 + yForward(:,2).^2 + yForward(:,3).^2 + yForward(:,4).^2 + 1e-8);
            
            y = yForward ./ max(1,cat(2,norm/lambda,norm/lambda,norm/lambda,norm/lambda));
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
    
    %clear comand output
    printToCmd( reverseStr,'');
    
    x = reshape(x,[size(image1),2]);
    y = reshape(y,[size(image1),4]);
        
end

