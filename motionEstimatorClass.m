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
        
        doGradientConstancy
        gradConstancyTermNumber
        
        dataTerm
        regularizerTerm

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
                obj.maxIt = 10000;
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
            
            if (exist('dataTerm','var'))
                obj.dataTerm = dataTerm;
            else
                obj.dataTerm = 'L1';
            end
            
            if (exist('doGradientConstancy','var'))
                obj.doGradientConstancy = doGradientConstancy;
            else
                obj.doGradientConstancy = 0;
            end
            
            if (exist('regularizerTerm','var'))
                obj.regularizerTerm = regularizerTerm;
            else
                obj.regularizerTerm = 'Huber';
            end
        end
        
        function init(obj)
            

            %a = diag(mask);
            %cut = find(a>1e-5, 1, 'first');
            %mask = mask(cut:100-cut+1, cut:100-cut+1);
            %mask = mask/sum(mask(:));
            
            %generate list of steps
            obj.steps = 1;
            i = 2;
            while i < Inf
                if (obj.steplength*obj.steps(i-1) * obj.dims{1}(1) > 25 && obj.steplength*obj.steps(i-1)*obj.dims{1}(2) > 25)
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
                % generate smoothing mask:
                smoothSigma = 1/sqrt(2*obj.steplength);
                mask = fspecial('gaussian', [99 99], smoothSigma);
                
                clear main;
                main = flexBox;
                main.params.maxIt = obj.maxIt;
                main.params.tryCPP = 1;
                
                tmpDims = obj.dims{i}(1:2);

                for j=1:obj.dims{i}(3)-1
                    %add two primal variables numbers are 2j-1 and 2j
                    numP1 = main.addPrimalVar(tmpDims);
                    numP2 = main.addPrimalVar(tmpDims);
                    
                    main.addTerm(emptyDataTerm(),numP1);
                    main.addTerm(emptyDataTerm(),numP2);
                    
                    if (i==1)
                        uTmp1{i,j} = obj.imageSequence(:,:,j);
                        uTmp2{i,j} = obj.imageSequence(:,:,j+1);
                    else
                        uTmp1{i,j} = imfilter(uTmp1{i-1,j}, mask, 'replicate');
                        uTmp2{i,j} = imfilter(uTmp2{i-1,j}, mask, 'replicate');
                    end

                    uTmp1{i,j} = imresize(uTmp1{i,j},tmpDims,'bicubic');
                    uTmp2{i,j} = imresize(uTmp2{i,j},tmpDims,'bicubic');

                    %add optical flow data term
                    if (strcmp(obj.dataTerm,'L2'))
                        main.addTerm(L2opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j}),[numP1,numP2]);
                        obj.flowTermNumber(j) = numel(main.duals);
                        
                        if (obj.doGradientConstancy)
                            main.addTerm(L2opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j},'termType','gradientConstancy','constancyDimension',1),[numP1,numP2]);
                            obj.gradConstancyTermNumber(j,1) = numel(main.duals);
                            main.addTerm(L2opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j},'termType','gradientConstancy','constancyDimension',2),[numP1,numP2]);
                            obj.gradConstancyTermNumber(j,2) = numel(main.duals);
                        end
                    else %default is L1
                        main.addTerm(L1opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j}),[numP1,numP2]);
                        obj.flowTermNumber(j) = numel(main.duals);
                        
                        if (obj.doGradientConstancy)
                            main.addTerm(L1opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j},'termType','gradientConstancy','constancyDimension',1),[numP1,numP2]);
                            obj.gradConstancyTermNumber(j,1) = numel(main.duals);
                            main.addTerm(L1opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j},'termType','gradientConstancy','constancyDimension',2),[numP1,numP2]);
                            obj.gradConstancyTermNumber(j,2) = numel(main.duals);
                        end
                    end
                    
                    %obj.flowTermNumber(j) = numel(main.duals);

                    %add regularizers - one for each component
                    if (strcmp(obj.regularizerTerm,'L2'))
                        main.addTerm(L2gradient(obj.alpha,tmpDims),numP1);
                        main.addTerm(L2gradient(obj.alpha,tmpDims),numP2);
                    elseif (strcmp(obj.regularizerTerm,'TV'))
                        main.addTerm(L1gradientIso(obj.alpha,tmpDims),numP1);
                        main.addTerm(L1gradientIso(obj.alpha,tmpDims),numP2);
                    else %default is Huber
                        main.addTerm(huberGradient(obj.alpha,tmpDims,0.01),numP1);
                        main.addTerm(huberGradient(obj.alpha,tmpDims,0.01),numP2);
                    end

                end
                
                obj.listFlexbox{i} = main;
            end
            
            if (obj.verbose > 0)
                disp('Initialization finished');
            end
        end
        
        function resetImages(obj,imageSequence)
            obj.imageSequence = imageSequence;
            
            % generate smoothing mask:
            smoothSigma = 1/sqrt(2*obj.steplength);
            mask = fspecial('gaussian', [99 99], smoothSigma);
            
            for i=1:numel(obj.steps)
                tmpDims = obj.dims{i}(1:2);
                
                for j=1:obj.dims{i}(3)-1
                    if (i==1)
                        uTmp1{i,j} = obj.imageSequence(:,:,j);
                        uTmp2{i,j} = obj.imageSequence(:,:,j+1);
                    else
                        uTmp1{i,j} = imfilter(uTmp1{i-1,j}, mask, 'replicate');
                        uTmp2{i,j} = imfilter(uTmp2{i-1,j}, mask, 'replicate');
                    end

                    uTmp1{i,j} = imresize(uTmp1{i,j},tmpDims,'bicubic');
                    uTmp2{i,j} = imresize(uTmp2{i,j},tmpDims,'bicubic');
                    
                    %add optical flow data term
                    if (strcmp(obj.dataTerm,'L2'))
                        obj.listFlexbox{i}.duals{obj.flowTermNumber(j)} = L2opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j});
                        
                        if (obj.doGradientConstancy)
                            obj.listFlexbox{i}.duals{obj.gradConstancyTermNumber(j,1)} = L2opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j},'termType','gradientConstancy','constancyDimension',1);
                            obj.listFlexbox{i}.duals{obj.gradConstancyTermNumber(j,2)} = L2opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j},'termType','gradientConstancy','constancyDimension',2);
                        end
                    else %default is L1
                        obj.listFlexbox{i}.duals{obj.flowTermNumber(j)} = L1opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j});
                        
                        if (obj.doGradientConstancy)
                            obj.listFlexbox{i}.duals{obj.gradConstancyTermNumber(j,1)} = L1opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j},'termType','gradientConstancy','constancyDimension',1);
                            obj.listFlexbox{i}.duals{obj.gradConstancyTermNumber(j,2)} = L1opticalFlowTerm(1,uTmp1{i,j},uTmp2{i,j},'termType','gradientConstancy','constancyDimension',2);
                        end
                    end
                end
            end
        end

        function v=getResult(obj)
            for j=1:obj.dims{1}(3)-1
                v(:,:,j,1) = obj.listFlexbox{1}.getPrimal(2*j-1);
                v(:,:,j,2) = obj.listFlexbox{1}.getPrimal(2*j);
            end
            %attach one empty slice
            v(:,:,obj.dims{1}(3),:) = 0;
        end
        
        function runLevel(obj,i)
            if (obj.verbose > 0)%on light verbose
                disp(['Running level #',num2str(i)])
            end
            
            for warps = 1:obj.numberOfWarps
                %warp data terms
                for j=1:obj.dims{i}(3)-1
                    obj.listFlexbox{i}.duals{obj.flowTermNumber(j)}.warpDataterm(obj.listFlexbox{i}.x{2*j-1},obj.listFlexbox{i}.x{2*j});
                    
                    if (obj.doGradientConstancy)
                        obj.listFlexbox{i}.duals{obj.gradConstancyTermNumber(j,1)}.warpDataterm(obj.listFlexbox{i}.x{2*j-1},obj.listFlexbox{i}.x{2*j});
                        obj.listFlexbox{i}.duals{obj.gradConstancyTermNumber(j,2)}.warpDataterm(obj.listFlexbox{i}.x{2*j-1},obj.listFlexbox{i}.x{2*j});
                    end
                end

                obj.listFlexbox{i}.runAlgorithm;
                
                if (obj.medianFiltering)
                    for j=1:obj.dims{i}(3)-1

                        obj.listFlexbox{i}.x{2*j-1} = medfilt2(obj.listFlexbox{i}.getPrimal(2*j-1), [5 5],'symmetric');
                        obj.listFlexbox{i}.x{2*j} = medfilt2(obj.listFlexbox{i}.getPrimal(2*j), [5 5],'symmetric');

                        obj.listFlexbox{i}.x{2*j-1} = obj.listFlexbox{i}.x{2*j-1}(:);
                        obj.listFlexbox{i}.x{2*j} = obj.listFlexbox{i}.x{2*j}(:);

                        if (obj.verbose > 1)%on massive verbose
                            clear flowField;
                            flowField(:,:,1) = obj.listFlexbox{i}.getPrimal(2*j-1);
                            flowField(:,:,2) = obj.listFlexbox{i}.getPrimal(2*j);

                            figure(300+j);imagesc(flowToColorV2(flowField));axis image;title(['Level #',num2str(i),' : Field #',num2str(j)]);drawnow;

                        end
                    end
                end
            end
            

        end
        
        function runPyramid(obj)
            if (obj.verbose > 0)%on light verbose
                disp(['Starting flow estimation pyramid'])
            end
            
            for i=numel(obj.steps) : -1 : 1
                %warp previous
                if (i < numel(obj.steps))
                    %upsample fields
                    tmpDimsOld = obj.dims{i+1}(1:2);
                    tmpDimsNew = obj.dims{i}(1:2);
                    
                    for j=1:numel(obj.listFlexbox{i}.x)
                        obj.listFlexbox{i}.x{j} = reshape(obj.listFlexbox{i+1}.x{j},tmpDimsOld);
                        obj.listFlexbox{i}.x{j} = imresize(obj.listFlexbox{i}.x{j},tmpDimsNew,'nearest')*tmpDimsNew(mod(1+j,2)+1)/tmpDimsOld(mod(1+j,2)+1);
                        obj.listFlexbox{i}.x{j} = obj.listFlexbox{i}.x{j}(:);
                    end
                    
                    for j=1:numel(obj.listFlexbox{i}.y)
                        obj.listFlexbox{i}.y{j} = reshape(obj.listFlexbox{i+1}.y{j},tmpDimsOld);
                        obj.listFlexbox{i}.y{j} = imresize(obj.listFlexbox{i}.y{j},tmpDimsNew,'nearest');
                        obj.listFlexbox{i}.y{j} = obj.listFlexbox{i}.y{j}(:);
                    end
                end
                obj.runLevel(i);
            end
            
            if (obj.verbose > 0)
                disp('Finished');
            end
        end
        
    end
    
end


