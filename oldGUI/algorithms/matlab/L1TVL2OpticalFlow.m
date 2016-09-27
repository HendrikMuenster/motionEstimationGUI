function [x,y] = L1TVL2OpticalFlow(u1,u2,tol,alpha0,alpha1,varargin)
    %insert is a list of images, output a list of velocity fields
    %returns primal and dual variable
    vararginParser;
    
    nPx = numel(u1);
    dimsU = size(u1);
    
    gradTV = fdgradient(dimsU,[1 1]);
    gradTV = gradTV.G;
    
    partialX = gradTV(1:nPx,:);
    partialY = gradTV(nPx+1:end,:);
    
    clear gradTV;
    
    partialXT = partialX';
    partialYT = partialY';

    tau=1/sqrt(8);
    sigma = tau;
    
    gradCentral = cdgradient(dimsU,[1 1]);
    gradCentral = gradCentral.G;
    
    Dy3Dc = gradCentral(1:nPx,:);
    Dx3Dc = gradCentral(nPx+1:2*nPx,:);

    ux = Dx3Dc*u1(:);
    uy = Dy3Dc*u1(:);
    ut = u2(:)-u1(:);
    
    %cleanup
    clear gradCentral Dy3Dc Dx3Dc;

    if (~exist('x','var'))
        x = zeros(nPx,6);
    end
    if (~exist('y','var'))
        y = zeros(nPx,4);
    end
    if (~exist('maxIt','var'))
        maxIt = 50000;
    end
    if (~exist('typeNorm','var'))
        %1 anisotropic
        %2 isotropic
        %3 fully isotropic
        typeNorm = 3;
    end
    
    Kty = x;
    Kx = y;

    aQuad = max(1e-9,ux.^2+uy.^2);
    tauAQuad = tau*aQuad;
    
    teiler1 = ux./aQuad;
    teiler2 = uy./aQuad;
    
    tauUx = tau*ux;
    tauUy = tau*uy;
    
    reverseStr = [];
    
    iterations = 1;err=1;
    while err>tol && iterations < maxIt
        xOld = x;
        yOld = y;
        
        if (mod(iterations,100)==1)
            if (exist('gtFlow','var'))
                v1 = reshape(x(:,1),dimsU);
                v2 = reshape(x(:,2),dimsU);
                
                GTerr = absoluteError( cat(3,v1,v2),gtFlow );
                
                reverseStr = printToCmd( reverseStr,sprintf('L1-TV-L2 Optical Flow, Iteration: #%d : Residual %f, AEE : %f',iterations,err,GTerr) );
            else
                reverseStr = printToCmd( reverseStr,sprintf('L1-TV-L2 Optical Flow, Iteration: #%d : Residual %f',iterations,err) );
            end
        end
        
        KtyOld = Kty;
        
        Kty(:,1) = partialXT*y(:,1) + partialYT*y(:,2);
        Kty(:,2) = partialXT*y(:,3) + partialYT*y(:,4);
        Kty(:,3) = -y(:,1);
        Kty(:,4) = -y(:,2);
        Kty(:,5) = -y(:,3);
        Kty(:,6) = -y(:,4);
        
        xHat = x - tau*Kty;

        rho = ut + ux.*xHat(:,1) + uy.*xHat(:,2);

        c1 = rho<= -tauAQuad;
        c2 = rho>= tauAQuad;
        c1c2 = c1-c2;
        c3 = 1-c1-c2;

        %primal prox
        x(:,1) = xHat(:,1) + tauUx.*c1c2 - teiler1 .* rho .* c3;
        x(:,2) = xHat(:,2) + tauUy.*c1c2 - teiler2 .* rho .* c3;
        x(:,3:6) = xHat(:,3:6) / (1+tau*alpha1);
        
        KxOld = Kx;
        
        Kx(:,1) = partialX*x(:,1) - x(:,3);
        Kx(:,2) = partialY*x(:,1) - x(:,4);
        Kx(:,3) = partialX*x(:,2) - x(:,5);
        Kx(:,4) = partialY*x(:,2) - x(:,6);

        %dual prox
        if (typeNorm == 1)
            y = max(-alpha0,min(alpha0,y + sigma*(2*Kx-KxOld)));
        elseif (typeNorm == 2)
            yTilde = y + sigma*(2*Kx-KxOld);
            
            norm1 = sqrt(yTilde(:,1).^2 + yTilde(:,2).^2);
            norm2 = sqrt(yTilde(:,3).^2 + yTilde(:,4).^2);
            
            y(:,1) = yTilde(:,1) ./ max(1,norm1/alpha0);
            y(:,2) = yTilde(:,2) ./ max(1,norm1/alpha0);
            y(:,3) = yTilde(:,3) ./ max(1,norm2/alpha0);
            y(:,4) = yTilde(:,4) ./ max(1,norm2/alpha0);
        elseif (typeNorm == 3)
            yTilde = y + sigma*(2*Kx-KxOld);
            
            norm1 = sqrt(yTilde(:,1).^2 + yTilde(:,2).^2 + yTilde(:,3).^2 + yTilde(:,4).^2);
            
            tmpDiv = max(1,norm1/alpha0);
            
            y(:,1) = yTilde(:,1) ./ tmpDiv;
            y(:,2) = yTilde(:,2) ./ tmpDiv;
            y(:,3) = yTilde(:,3) ./ tmpDiv;
            y(:,4) = yTilde(:,4) ./ tmpDiv;
        end
        
        xDiff = xOld - x;
        yDiff = yOld - y;

        %primal residual
        p = xDiff/tau - (KtyOld-Kty);
        
        %dual residual
        d = yDiff/sigma - (KxOld-Kx);

        p=sum(abs(p(:)));
        d=sum(abs(d(:)));
        
        err = (p+d) / (2*nPx);

        iterations = iterations + 1;
    end
    
    %v = cat(3,reshape(x(:,1),dimsU),reshape(x(:,2),dimsU));
    
    %clear comand output
    printToCmd( reverseStr,'');
        
end

