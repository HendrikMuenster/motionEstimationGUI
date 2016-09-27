function [x,y] = motionEstimationPyramidFlexBox(u,dimsU,tol,alpha,algorithmName,numDualVars,varargin)
    nPx = prod(dimsU);

    if (numel(varargin) == 1)
        varargin = [varargin{:}];
    end
    vararginParser;

	initVar('verbose',0);
    initVar('numPrimalVars',2);
    initVar('x',zeros([dimsU,numPrimalVars]));
    initVar('y',zeros([dimsU,numDualVars]));
    initVar('maxIt',10000);
    initVar('typeNorm',4);
    initVar('stepsize',[1 1]);
    initVar('discretization',1);
    initVar('numsteps',200);
    initVar('steplength',0.8);
    initVar('numberOfWarps',3);
    initVar('huberEpsilon',0.01);
    
    
    smoothSigma = 1/sqrt(2*steplength);
    initVar('doGaussianSmoothing',1);
    initVar('medianFiltering',1);
    initVar('medianFilteringTerm',0);
    initVar('adjustStepsize',1);
    initVar('doGradientConstancy',1);
    
    %option to define useCPP from outside
    if (~exist('useCPP','var'))
        if (exist([algorithmName,'CPP'],'file'))
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
        
        checkImage = imresize(u(:,:,1),steps(i),'bicubic');

        newHeight = size(checkImage,1);
        newWidth = size(checkImage,2);

        newSize = [newHeight newWidth];
        
        uTmp = [];
        xTmp = [];
        yTmp = [];
        
        for j=1:dimsU(3)
            if (doGaussianSmoothing)
                uTmp2 = u(:,:,j);
                uTmp2 = imfilter(uTmp2, mask, 'replicate');

                uTmp(:,:,j) = imresize(uTmp2,newSize,'bicubic');
            else
                uTmp(:,:,j) = imresize(u(:,:,j),newSize,'bicubic');
            end
            
            if (j < dimsU(3))
                %rescale velocity field
                xTmp(:,:,j,1) = imresize(x(:,:,j,1),newSize,'nearest')*newHeight/oldHeight;
                xTmp(:,:,j,2) = imresize(x(:,:,j,2),newSize,'nearest')*newWidth/oldWidth;
                
                for k=1:numDualVars
                    yTmp(:,:,j,k) = imresize(y(:,:,j,k),newSize,'nearest');
                end
            end
        end
        
        if (adjustStepsize)
            %adjust steplength to reference grid of max length 256
            maxDim = max(dimsU(1),dimsU(2));
            
            stepsize(1) = 512/maxDim * dimsU(1) / newSize(1);
            stepsize(2) = 512/maxDim * dimsU(2) / newSize(2);
        end
        
        %init flexbox for this level
        
        main = flexBox;
        main.params.maxIt = maxIt;
        main.params.tryCPP = 1;
        
        for j=1:dimsU(3)-1
            %add two primal variables numbers are 2j-1 and 2j
            main.addPrimalVar(newSize);
            main.addPrimalVar(newSize);
            
            %add optical flow data term
            flowTermNumber(j) = main.addTerm(L1opticalFlowTerm(1,uTmp(:,:,j),uTmp(:,:,j+1)),[2*j-1,2*j]);
            
            if (doGradientConstancy)
                gradConstancyTermNumber1(j) = main.addTerm(L1opticalFlowTerm(1,uTmp(:,:,j),uTmp(:,:,j+1),'termType','gradientConstancy','constancyDimension',1),[2*j-1,2*j]);
                gradConstancyTermNumber2(j) = main.addTerm(L1opticalFlowTerm(1,uTmp(:,:,j),uTmp(:,:,j+1),'termType','gradientConstancy','constancyDimension',2),[2*j-1,2*j]);
            end
            
            %add regularizers - one for each component
            
            %main.addTerm(L1gradientIso(alpha,newSize,'stepsize',stepsize),2*j-1);
            %main.addTerm(L1gradientIso(alpha,newSize,'stepsize',stepsize),2*j);
            main.addTerm(huberGradient(alpha,newSize,huberEpsilon,'stepsize',stepsize),2*j-1);
            main.addTerm(huberGradient(alpha,newSize,huberEpsilon,'stepsize',stepsize),2*j);
            
            if (medianFilteringTerm)
                [li,lj,lv] = createMedianFilter(size(uTmp(:,:,j)),5);

                operator = sparse(li,lj,lv);
                main.addTerm(L1operatorAniso(0.01,1,operator),2*j-1);
                main.addTerm(L1operatorAniso(0.01,1,operator),2*j);
            end
        end
        
        if (i>1)
            main.x{1} = reshape(xTmp(:,:,j,1),[numel(xTmp(:,:,j,1)),1]);
            main.x{2} = reshape(xTmp(:,:,j,2),[numel(xTmp(:,:,j,2)),1]);
            for k=1:numDualVars
                main.y{numDualVars * (j-1) + k} = reshape(yTmp(:,:,j,k),[numel(yTmp(:,:,j,k)),1]);
            end
        end

        for warps = 1:numberOfWarps
            
            %warp data terms
            for j=1:dimsU(3)-1
                main.duals{flowTermNumber(j)}.warpDataterm(xTmp(:,:,j,1),xTmp(:,:,j,2));
                if (doGradientConstancy)
                    main.duals{gradConstancyTermNumber1(j)}.warpDataterm(xTmp(:,:,j,1),xTmp(:,:,j,2));
                    main.duals{gradConstancyTermNumber2(j)}.warpDataterm(xTmp(:,:,j,1),xTmp(:,:,j,2));
                end
            end
            
            main.runAlgorithm;
            
            %do median filtering
            for j=1:dimsU(3)-1
                xTmp(:,:,j,1) = main.getPrimal(2*j-1);
                xTmp(:,:,j,2) = main.getPrimal(2*j);
                 
                if (medianFiltering)
                    xTmp(:,:,j,1) = medfilt2(xTmp(:,:,j,1), [5 5],'symmetric');
                    xTmp(:,:,j,2) = medfilt2(xTmp(:,:,j,2), [5 5],'symmetric');
                end
                
                for k=1:numDualVars
                    yTmp(:,:,j,k) = reshape(main.getDual( numDualVars * (j-1) + k),size(xTmp(:,:,j,1)));
                end
            end
            
            flowField = squeeze(xTmp(:,:,1,:));axis image;drawnow;
        

            %figure(4);imagesc(flowToColorV2(flowField));axis image;drawnow;
            
            
        end
        
        %get results
        x = xTmp;
        y = yTmp;
        
        flowField = squeeze(x(:,:,1,:));
        

        figure(4);imagesc(flowToColorV2(flowField));axis image;drawnow;
        
    end
    
    %add zero frames
    x(:,:,end+1,:) = 0;
    y(:,:,end+1,:) = 0;
   
    printToCmd(strInner,'');
end

