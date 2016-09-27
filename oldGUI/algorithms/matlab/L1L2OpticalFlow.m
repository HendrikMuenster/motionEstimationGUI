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
function [x,y] = L1L2OpticalFlow(image1,image2,tol,lambda,varargin)
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
    tau = 1/sqrt(8);
    sigma = tau;
    
    betaQuad = max(1e-9,ux.^2+uy.^2);
    tauBetaQuad = tau*betaQuad;
    
    teiler(:,1) = ux ./ betaQuad;
    teiler(:,2) = uy ./ betaQuad;
    
    tauGrad(:,1) = tau*ux;
    tauGrad(:,2) = tau*uy;
    
    if (~exist('x','var'))
        x = zeros(nPx,2);
    end
    if (~exist('y','var'))
        y = zeros(nPx,4);
    end
    if (~exist('maxIt','var'))
        maxIt = 50000;
    end
    
    Kx = y;
    Kty = x;

    reverseStr = [];
    
    iterations = 1;err=1;
    while err>tol && iterations < maxIt
        yOld = y;
        xOld = x;
        
        if (mod(iterations,100)==1)
            reverseStr = printToCmd( reverseStr,sprintf('L2-L2 Optical Flow, Iteration: #%d : Residual %f',iterations,err) );
            
            %xR = reshape(x,[size(image1),2]);
            %figure(4);clf;imagesc(colourfulOrientationPlot(xR(:,:,1),xR(:,:,2),5));drawnow
        end
        
        KtyOld = Kty;
        
        Kty(:,1) = Dx2D'*y(:,1) + Dy2D'*y(:,2);
        Kty(:,2) = Dx2D'*y(:,3) + Dy2D'*y(:,4);
        
        xHat = x - tau*Kty;
        
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
        y = lambda/(sigma+lambda)* (y + sigma*(2*Kx-KxOld)) ;

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

