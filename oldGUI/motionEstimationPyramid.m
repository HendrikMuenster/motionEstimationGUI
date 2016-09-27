function [x,y] = motionEstimationPyramid(u,dimsU,tol,alpha,algorithmName,numDualVars,varargin)
    nPx = prod(dimsU);

    if (numel(varargin) == 1)
        varargin = [varargin{:}];
    end
    vararginParser;
    
    initVar('useCPP',0);
    initVar('useCudaCPP',0);

	initVar('verbose',0);
    initVar('numPrimalVars',2);
    initVar('x',zeros([dimsU,numPrimalVars]));
    initVar('y',zeros([dimsU,numDualVars]));
    initVar('maxIt',10000);
    initVar('typeNorm',4);
    initVar('huberEpsilon',0.01);
    
    initVar('stepsize',[1 1 1]);
    initVar('discretization',1);
    initVar('numsteps',200);
    initVar('steplength',0.8);
    initVar('numberOfWarps',3);
    
    
    
    smoothSigma = 1/sqrt(2*steplength);
    initVar('doGaussianSmoothing',1);
    initVar('medianFiltering',1);
    initVar('adjustStepsize',0);
    
    initVar('gradientConstancy',0);
    
     %option to define useCPP from outside
    if (useCPP ~= 2)
        if (exist([algorithmName,'CudaCPP'],'file'))
            useCudaCPP = 1;
        elseif (exist([algorithmName,'CPP'],'file'))
            useCPP = 1;
        else
            useCPP = 0;
        end
    end
    
    
    % generate smoothing mask:
    mask = fspecial('gaussian', [100, 100], smoothSigma);
    a = diag(mask);
    cut = find(a>1e-5, 1, 'first');
    mask = mask(cut:100-cut+1, cut:100-cut+1);
    mask = mask/sum(mask(:));

    %automatic step length generation
    steps = 1;
    for i=2:numsteps
        if (steplength*steps(i-1)*size(u,1) > 10 && steplength*steps(i-1)*size(u,2) > 10)
            steps(i) = steplength*steps(i-1);
        else
            break;
        end
    end
    steps = fliplr(steps);
    
    
    strInner = [];

    for i=1:numel(steps)
        if (verbose > 0)
            disp(['Starting zoom-level ',num2str(steps(i))]);
        end

        %no variables set
        if (i == 1)
            oldHeight = size(u,1);
            oldWidth = size(u,2);
        else
            oldHeight = size(uTmp,1);
            oldWidth = size(uTmp,2);
        end
        
        checkImage = imresize(u(:,:,1),steps(i),'cubic');

        newHeight = size(checkImage,1);
        newWidth = size(checkImage,2);

        newSize = [newHeight newWidth];

        uTmp = [];
        xTmp = [];
        yTmp = [];

        %rescale all variable to current level
        for j=1:dimsU(3)
            if (doGaussianSmoothing)
                uTmp2 = u(:,:,j);
                uTmp2 = imfilter(uTmp2, mask, 'replicate');

                uTmp(:,:,j) = imresize(uTmp2,newSize,'cubic');
            else
                uTmp(:,:,j) = imresize(u(:,:,j),newSize,'cubic');
            end
            
            if (j < dimsU(3))
                %rescale velocity field
                xTmp(:,:,j,1) = imresize(x(:,:,j,1),newSize,'nearest')*newHeight/oldHeight;
                xTmp(:,:,j,2) = imresize(x(:,:,j,2),newSize,'nearest')*newWidth/oldWidth;
                
                for k=3:numPrimalVars
                    xTmp(:,:,j,k) = imresize(x(:,:,j,k),newSize,'nearest');
                end

                for k=1:numDualVars
                    yTmp(:,:,j,k) = imresize(y(:,:,j,k),newSize,'nearest');
                end
            end
        end
        
        x = [];
        y = [];
        
        if (adjustStepsize)
            %adjust steplength to reference grid of max length 256
            maxDim = max(dimsU(1),dimsU(2));
            
            stepsize(2) = 512/maxDim * dimsU(1) / newSize(1);
            stepsize(3) = 512/maxDim * dimsU(2) / newSize(2);
            stepsize(2) = dimsU(1) / newSize(1);
            stepsize(3) = dimsU(2) / newSize(2);
        end
        %stepsize*newSize(1)

        for j=1:dimsU(3)-1

            if (useCPP || useCudaCPP)
                if (useCudaCPP)
                    cudaString = 'Cuda';
                else
                    cudaString = '';
                end
                
                xT = xTmp(:,:,j,:);
                yT = yTmp(:,:,j,:);
            
                if (exist('alpha1','var'))
                    [xRes,yRes] = eval([algorithmName,'CPP','(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,alpha1,maxIt,typeNorm,xT,yT,stepsize,discretization,numberOfWarps)']);
                    
                    xRes = reshape(xRes,[newHeight,newWidth,numPrimalVars]);
                    for k=1:numPrimalVars
                        x(:,:,j,k) = xRes(:,:,k);
                    end
                else
                    if (strcmp(algorithmName,'L2TVBregOpticalFlow'))
                        [x(:,:,j,1),x(:,:,j,2),yRes] = eval(['L2TVBregOpticalFlowCPP','(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,maxIt,numBreg)']);
                    elseif (strcmp(algorithmName,'L1TVOpticalFlowNonlinear'))
                        [x(:,:,j,1),x(:,:,j,2),yRes] = eval([algorithmName,cudaString,'CPP','(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,maxIt,typeNorm,xT,yT,stepsize,discretization,numberOfWarps,huberEpsilon,gradientConstancy)']);
                    else
                        [x(:,:,j,1),x(:,:,j,2),yRes] = eval([algorithmName,cudaString,'CPP','(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,maxIt,typeNorm,xT,yT,stepsize,discretization,numberOfWarps)']);
                    end
                    
                    
                end

                yRes = reshape(yRes,[newHeight,newWidth,numDualVars]);
                for k=1:numDualVars
                    y(:,:,j,k) = yRes(:,:,k);
                end
            else
                clear xT yT;
                
                for k=1:numPrimalVars
                    xT(:,k) = reshape(xTmp(:,:,j,k),[prod(newSize),1]);
                end
                for k=1:numDualVars
                    yT(:,k) = reshape(yTmp(:,:,j,k),[prod(newSize),1]);
                end
                
                if (exist('alpha1','var'))
                    [xRes,yRes] = eval([algorithmName,'(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,alpha1,''x'',xT,''y'',yT,''maxIt'',maxIt)']);
                else
                    [xRes,yRes] = eval([algorithmName,'(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,''x'',xT,''y'',yT,''maxIt'',maxIt,''stepsize'',stepsize,''huberEpsilon'',huberEpsilon,''typeNorm'',typeNorm)']);
                end
                
                for k=1:numPrimalVars
                    x(:,:,j,k) = xRes(:,:,k);
                end
                for k=1:numDualVars
                    y(:,:,j,k) = yRes(:,:,k);
                end
            end
            
            if (medianFiltering)
                x(:,:,j,1) = medfilt2(x(:,:,j,1), [5 5],'symmetric');
                x(:,:,j,2) = medfilt2(x(:,:,j,2), [5 5],'symmetric');
            end

        end
        
        flowField = squeeze(x(:,:,1,:));
        

        figure(4);imagesc(flowToColorV2(flowField));axis image;drawnow;
        %pause
    end
    
    %add zero frames
    x(:,:,end+1,:) = 0;
    y(:,:,end+1,:) = 0;
   
    printToCmd(strInner,'');
end

