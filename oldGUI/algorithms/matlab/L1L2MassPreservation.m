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
function [x,y] = L1L2MassPreservation(image1,image2,tol,lambda,varargin)
    %insert is a list of images, output a list of velocity fields
    %returns primal and dual variable
    vararginParser;
    
    nPx = numel(image1);
    
    [ Dy2Dc,Dx2Dc ] = generateCentralGradient2D( size(image1),[1 1]); 
    %switch Dx and Dy
    
    %Dx*u
    Dx2Dc = Dx2Dc * sparse(1:nPx,1:nPx,image1(:));
    %Dy*u
    Dy2Dc = Dy2Dc * sparse(1:nPx,1:nPx,image1(:));
    
    ut = image2(:)-image1(:);
    
    [ Dx2D,Dy2D ] = generateForwardGradient2D( size(image1),[1 1]);
    
    nPx = numel(image1);
    tau = 1/sqrt(8);
    sigma = tau;
    
    sigUt = sigma*ut;
    
    if (~exist('x','var'))
        x = zeros(nPx,2);
    end
    if (~exist('y','var'))
        y = zeros(nPx,5);
    end
    if (~exist('maxIt','var'))
        maxIt = 50000;
    end
    if (~exist('typeNorm','var'))
        typeNorm = 3;
    end
    
    Kx = y;
    Kty = x;
    
    reverseStr = [];
    
    iterations = 1;err=1;
    while err>tol && iterations < maxIt
        yOld = y;
        xOld = x;
        
        if (mod(iterations,100)==1)
            reverseStr = printToCmd( reverseStr,sprintf('L2-L2 Mass Preservation, Iteration: #%d : Residual %f',iterations,err) );
            
            %xR = reshape(x,[size(image1),2]);
            %figure(1);clf;imagesc(colourfulOrientationPlot(xR(:,:,1),xR(:,:,2),5));drawnow
        end
        
        KtyOld = Kty;%K(y_{k})
        
        Kty(:,1) = Dx2D'*y(:,1) + Dy2D'*y(:,2) + Dx2Dc'*y(:,5);
        Kty(:,2) = Dx2D'*y(:,3) + Dy2D'*y(:,4) + Dy2Dc'*y(:,5);
        
        %primal prox
        x = x - tau*Kty;
        
        KxOld = Kx;
        Kx(:,1) = Dx2D*x(:,1);
        Kx(:,2) = Dy2D*x(:,1);
        Kx(:,3) = Dx2D*x(:,2);
        Kx(:,4) = Dy2D*x(:,2);
        Kx(:,5) = Dx2Dc*x(:,1) + Dy2Dc*x(:,2);
        
        yForward = y + sigma*(2*Kx-KxOld);
            
        %dual prox
        y(:,1:4) = lambda/(sigma+lambda)* (y(:,1:4) + sigma*(2*Kx(:,1:4)-KxOld(:,1:4))) ;
        
        %mass preservation part
        y(:,5) = max(-1,min(1,yForward(:,5) + sigUt));
        
        xDiff = x - xOld;
        yDiff = y - yOld;

        %primal residul
        p = xDiff/tau - (KtyOld-Kty);
        %dual residual
        d = yDiff/sigma - (KxOld-Kx);

        p=sum(abs(p(:)));
        d=sum(abs(d(:)));

        err = (p+d) / nPx;

        iterations = iterations + 1;
    end
    
    %clear comand output
    printToCmd(reverseStr,'');
    
    x = reshape(x,[size(image1),2]);
    y = reshape(y,[size(image1),5]);
end

