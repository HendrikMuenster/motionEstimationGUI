classdef motionEstimatorClass < handle

    properties
        verbose
        tol
        alpha
        maxIt
        steplength
        steps
        numberOfWarps
        doGaussianSmoothing
        medianFiltering
        imageSequence
        dims
        flowTermNumber

        listFlexbox
    end
    
    methods
        function obj = motionEstimatorClass(imageSequence,tol,alpha,varargin)
            vararginParser;
            
            obj.imageSequence = imageSequence;
            obj.alpha = alpha;
            obj.tol = tol;
            obj.dims{1} = size(imageSequence);
            
            if (exist('verbose','var'))
                obj.verbose = verbose;
            else
                obj.verbose = 0;
            end
            
            if (exist('steplength','var'))
                obj.steplength = steplength;
            else
                obj.steplength = 0.8;
            end
            
            if (exist('maxIt','var'))
                obj.maxIt = maxIt;
            else
                obj.maxIt = 100;
            end
            
            if (exist('numberOfWarps','var'))
                obj.numberOfWarps = numberOfWarps;
            else
                obj.numberOfWarps = 3;
            end
            
            if (exist('doGaussianSmoothing','var'))
                obj.doGaussianSmoothing = doGaussianSmoothing;
            else
                obj.doGaussianSmoothing = 1;
            end
            
            if (exist('medianFiltering','var'))
                obj.medianFiltering = medianFiltering;
            else
                obj.medianFiltering = 1;
            end
        end
        
        function init(obj)
            
            % generate smoothing mask:
            smoothSigma = 1/sqrt(2*obj.steplength);
            mask = fspecial('gaussian', [100, 100], smoothSigma);
            a = diag(mask);
            cut = find(a>1e-5, 1, 'first');
            mask = mask(cut:100-cut+1, cut:100-cut+1);
            mask = mask/sum(mask(:));
            
            %generate list of steps
            obj.steps = 1;
            i = 2;
            while i < Inf
                if (obj.steplength*obj.steps(i-1) * obj.dims{1}(1) > 10 && obj.steplength*obj.steps(i-1)*obj.dims{1}(2) > 10)
                    obj.steps(i) = obj.steplength*obj.steps(i-1);
                    
                    obj.dims{i} = round(obj.dims{1} * obj.steps(i));
                    obj.dims{i}(3) = obj.dims{1}(3);
                    
                    i = i + 1;
                else
                    break;
                end
            end
            
            %% create flexBox obejects
            
            for i=1:numel(obj.steps)
                main = flexBox;
                main.params.maxIt = obj.maxIt;
                main.params.tryCPP = 1;
                
                tmpDims = obj.dims{i}(1:2);
                
%                 for j=1:obj.dims{i}(3)
%                     figure(12+j);imagesc(obj.imageSequence(:,:,j));axis image;
%                     pause
%                 end

                for j=1:obj.dims{i}(3)-1
                    %add two primal variables numbers are 2j-1 and 2j
                    numP1 = main.addPrimalVar(tmpDims);
                    numP2 = main.addPrimalVar(tmpDims);
                    
                    main.addTerm(emptyDataTerm(),numP1);
                    main.addTerm(emptyDataTerm(),numP2);
                    
                    uTmp1 = imfilter(obj.imageSequence(:,:,j), mask, 'replicate');
                    uTmp2 = imfilter(obj.imageSequence(:,:,j+1), mask, 'replicate');

                    uTmp1 = imresize(uTmp1,tmpDims,'bicubic');
                    uTmp2 = imresize(uTmp2,tmpDims,'bicubic');

                    %add optical flow data term
                    main.addTerm(L1opticalFlowTerm(1,uTmp1,uTmp2),[numP1,numP2]);
                    obj.flowTermNumber(j) = numel(main.duals);

                    %add regularizers - one for each component
                    main.addTerm(huberGradient(obj.alpha,tmpDims,0.01),numP1);
                    main.addTerm(huberGradient(obj.alpha,tmpDims,0.01),numP2);
                    
                    %figure(12);imagesc(uTmp1);axis image;
                    %figure(13);imagesc(uTmp2);axis image;
                end
                
                obj.listFlexbox{i} = main;
            end
            
            
            disp('Initialization finished');
        end
        
        function resetImages(obj,imageSequence)
            obj.imageSequence = imageSequence;
            
            % generate smoothing mask:
            smoothSigma = 1/sqrt(2*obj.steplength);
            mask = fspecial('gaussian', [100, 100], smoothSigma);
            a = diag(mask);
            cut = find(a>1e-5, 1, 'first');
            mask = mask(cut:100-cut+1, cut:100-cut+1);
            mask = mask/sum(mask(:));
            
            for i=1:numel(obj.steps)
                
                tmpDims = obj.dims{i}(1:2);
                
                for j=1:obj.dims{1}(3)-1
                    uTmp1 = imfilter(obj.imageSequence(:,:,j), mask, 'replicate');
                    uTmp2 = imfilter(obj.imageSequence(:,:,j+1), mask, 'replicate');

                    uTmp1 = imresize(uTmp1,tmpDims,'bicubic');
                    uTmp2 = imresize(uTmp2,tmpDims,'bicubic');
                    
                    obj.listFlexbox{i}.duals{3*j-2} = L1opticalFlowTerm(1,uTmp1,uTmp2);
                end
            end
        end

        function v=getResult(obj)
            for j=1:obj.dims{1}(3)-1
                v(:,:,j,1) = obj.listFlexbox{1}.getPrimal(2*j-1);
                v(:,:,j,2) = obj.listFlexbox{1}.getPrimal(2*j);
            end
        end
        
        function runLevel(obj,i)
            for warps = 1:obj.numberOfWarps
                %warp data terms
                for j=1:obj.dims{i}(3)-1
                    obj.listFlexbox{i}.duals{obj.flowTermNumber(j)}.warpDataterm(obj.listFlexbox{i}.x{2*j-1},obj.listFlexbox{i}.x{2*j});
                end

                obj.listFlexbox{i}.runAlgorithm;

                if (obj.medianFiltering)
                    for j=1:obj.dims{i}(3)-1
                        obj.listFlexbox{i}.x{2*j-1} = medfilt2(obj.listFlexbox{i}.x{2*j-1}, [5 5],'symmetric');
                        obj.listFlexbox{i}.x{2*j} = medfilt2(obj.listFlexbox{i}.x{2*j}, [5 5],'symmetric');
                    
                    
                        clear flowField;
                        flowField(:,:,1) = obj.listFlexbox{i}.getPrimal(2*j-1);
                        flowField(:,:,2) = obj.listFlexbox{i}.getPrimal(2*j);

                        if (obj.verbose > 0)
                            figure(4+j);imagesc(flowToColorV2(flowField));axis image;drawnow;
                            
                        end
                    end
                end


            end
        end
        
        function runPyramid(obj)
            for i=numel(obj.steps) : -1 : 1
                %warp previous
                if (i < numel(obj.steps))
                    %upsample fields
                    tmpDimsOld = obj.dims{i+1}(1:2);
                    tmpDimsNew = obj.dims{i}(1:2);
                    
                    for j=1:numel(obj.listFlexbox{i}.x)
                        obj.listFlexbox{i}.x{j} = reshape(obj.listFlexbox{i+1}.x{j},tmpDimsOld);
                        obj.listFlexbox{i}.x{j} = imresize(obj.listFlexbox{i}.x{j},tmpDimsNew,'nearest')*tmpDimsNew(mod(1+j,2)+1)/tmpDimsOld(mod(1+j,2)+1);
                    end
                    
                    for j=1:numel(obj.listFlexbox{i}.y)
                        obj.listFlexbox{i}.y{j} = reshape(obj.listFlexbox{i+1}.y{j},tmpDimsOld);
                        obj.listFlexbox{i}.y{j} = imresize(obj.listFlexbox{i}.y{j},tmpDimsNew,'nearest');
                    end
                end
                
                obj.runLevel(i);
            end
            
            disp('Finished');
        end
        
    end
    
end


