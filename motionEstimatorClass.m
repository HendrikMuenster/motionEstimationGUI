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
        doWarping
        doGaussianSmoothing
        medianFiltering
        imageSequence
        dims
        
        flowTermNumber
        
        imageDiscretization
        
        doGradientConstancy
        gradConstancyTermNumber
        regularizerTermNumber
        
        dataTerm
        regularizerTerm

        flexboxObj
        
        %changes for 3D data
        is3D
        numberSpatialDims
        numberImages
        listImages
    end
    
    methods
        function obj = motionEstimatorClass(imageSequence,tol,alpha,varargin)
            vararginParser;
            
            obj.imageSequence = imageSequence;
            obj.alpha = alpha;
            obj.tol = tol;
            obj.dims{1} = size(imageSequence);
            
            if (ndims(imageSequence) == 4)
                obj.is3D = 1;
            else
                obj.is3D = 0;
            end
            
            obj.numberSpatialDims = ndims(obj.imageSequence) - 1;
            obj.numberImages = size(imageSequence);
            obj.numberImages = obj.numberImages(end);
            
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
            
            if (exist('doWarping','var'))
                obj.doWarping = doWarping;
                if (doWarping == 0)
                    obj.numberOfWarps = 1;
                end
            else
                obj.doWarping = 1;
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
            
            %discretization for the image derivatives.
            if (exist('imageDiscretization','var'))
                obj.imageDiscretization = imageDiscretization;
            else
                obj.imageDiscretization = 'interpolated';
            end
            
        end
        
        function init(obj)
            %% generate list of steps
            obj.steps = 1;
            i = 2;
            while i < Inf
                if (sum(obj.dims{i-1}(1:end-1) < 15) > 0)
                    break;
                else
                    obj.steps(i) = obj.steplength*obj.steps(i-1);
                    
                    obj.dims{i} = round(obj.dims{1} * obj.steps(i));

                    obj.dims{i}(end) = obj.dims{1}(end);
                    
                    i = i + 1;
                end
            end
            
            %% generate smoothing mask:
            %smoothSigma = 1/sqrt(2*obj.steplength);
            %mask = fspecial('gaussian', [99 99], smoothSigma);
            
            smooth_sigma      = sqrt(1/obj.steplength)/sqrt(2);
            mask = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);
            
            %% generate list of images
            for i=1:numel(obj.steps)
                for j=1:obj.numberImages
                    if (i==1)
                        idx = repmat({':'},1,obj.numberSpatialDims);idx{end+1} = j;
                        obj.listImages{i,j} = obj.imageSequence(idx{:});
                    else
                        tmpDimsOld = size(obj.listImages{i-1,j});
                        tmpDimsNew = obj.dims{i}(1:end-1);
                        
                        %create idx
                        for k=1:obj.numberSpatialDims
                            idxOld{k} = (0:(tmpDimsNew(k)-1) / (tmpDimsOld(k)-1) :tmpDimsNew(k)-1) + 1;
                            idxNew{k} = 1:tmpDimsNew(k);
                        end
                        
                        gridOld = cell(1,obj.numberSpatialDims);
                        gridNew = cell(1,obj.numberSpatialDims);
                        
                        [gridOld{:}] = ndgrid(idxOld{:});
                        [gridNew{:}] = ndgrid(idxNew{:});
                        
                        stringResize = '';
                        stringResize2 = '';
                        for k=1:numel(tmpDimsOld)
                            stringResize = [stringResize,',gridOld{',num2str(k),'}'];
                            stringResize2 = [stringResize2,',gridNew{',num2str(k),'}'];
                        end
                        stringResize = stringResize(2:end);
                        
                        obj.listImages{i,j} = imfilter(obj.listImages{i-1,j}, mask, 'replicate');
                        obj.listImages{i,j} = eval(['interpn(',stringResize,',obj.listImages{',num2str(i),',',num2str(j),'}',stringResize2,',''bicubic'');']);
                        
                        check1 = isnan(obj.listImages{i,j});
                        check2 = isinf(obj.listImages{i,j});

                        if ((sum(check1(:)) + sum(check2(:))) > 0)
                            error('Isinf + Isnan!!!')
                        end
                        
                    end
                end
            end

            
            %% create flexBox obeject
            obj.flexboxObj = flexBox;
            obj.flexboxObj.params.verbose = obj.verbose;
            obj.flexboxObj.params.maxIt = obj.maxIt;
            obj.flexboxObj.params.tryCPP = 1;
            
            %create initial flexBox object on coarsest scale
            tmpDims = obj.dims{1}(1:end-1);

            for j=1:obj.dims{1}(end) - 1
                %add primal variables
                for k=1:obj.numberSpatialDims
                    numP(k) = obj.flexboxObj.addPrimalVar(tmpDims);
                    %obj.flexboxObj.addTerm(emptyDataTerm(),numP(k));
                end
                
                if (strcmp(obj.dataTerm,'L2'))
                    flowTermObject = L2opticalFlowTerm(1,obj.listImages{1,j},obj.listImages{1,j+1},'discretization',obj.imageDiscretization);
                else %default is L1
                    flowTermObject = L1opticalFlowTerm(1,obj.listImages{1,j},obj.listImages{1,j+1},'discretization',obj.imageDiscretization);
                end
                
                obj.flexboxObj.addTerm(flowTermObject,numP);
                obj.flowTermNumber(j) = numel(obj.flexboxObj.duals);

                if (obj.doGradientConstancy)
                    for k=1:obj.numberSpatialDims
                        if (strcmp(obj.dataTerm,'L2'))
                            gradientConstancyTermObject = L2opticalFlowTerm(1,obj.listImages{1,j},obj.listImages{1,j+1},'termType','gradientConstancy','constancyDimension',k);
                        else %default is L1
                            gradientConstancyTermObject = L1opticalFlowTerm(1,obj.listImages{1,j},obj.listImages{1,j+1},'termType','gradientConstancy','constancyDimension',k);
                        end
                        
                        obj.flexboxObj.addTerm(gradientConstancyTermObject,numP);
                        obj.gradConstancyTermNumber(j,k) = numel(obj.flexboxObj.duals);
                    end
                end

                %create regularizer term
                if (strcmp(obj.regularizerTerm,'L2'))
                    regularizerObject = L2gradient(obj.alpha,tmpDims);
                elseif (strcmp(obj.regularizerTerm,'TV'))
                    regularizerObject = L1gradientIso(obj.alpha,tmpDims);
                else %default is Huber
                    regularizerObject = huberGradient(obj.alpha,tmpDims,0.01);
                end
                
                for k=1:numel(tmpDims)
                    obj.flexboxObj.addTerm(regularizerObject,numP(k));
                    obj.regularizerTermNumber(j,k) = numel(obj.flexboxObj.duals);
                end
            end

            if (obj.verbose > 0)
                disp('Initialization finished');
            end
        end
        
        function resetImages(obj,imageSequence)
            obj.imageSequence = imageSequence;
        end

        function v=getResult(obj)
            for j=1:obj.numberImages - 1
                for k=1:obj.numberSpatialDims
                    idx = repmat({':'},1,obj.numberSpatialDims);
                    idx{end+1} = j;
                    idx{end+1} = k;
                    
                    v(idx{:}) = obj.flexboxObj.getPrimal(obj.numberSpatialDims*(j-1) + k);
                end
            end
            %attach one empty slice
            idx = repmat({':'},1,obj.numberSpatialDims + 2);
            idx{end-1} = obj.numberImages;
            
            v(idx{:}) = 0;
        end
        
        function runLevel(obj,level)
            if (obj.verbose > 0)%on light verbose
                disp(['Running level #',num2str(level)])
            end
            
            %
            % adjust flexBox object for current level
            %
            
            %read old dimension
            tmpDimsOld = obj.flexboxObj.dims{1};
            tmpDimsNew = obj.dims{level}(1:end-1);
            
            %create idx
            for i=1:obj.numberSpatialDims
                idxOld{i} = (0:(tmpDimsNew(i)-1) / (tmpDimsOld(i)-1) :tmpDimsNew(i)-1) + 1;
                idxNew{i} = 1:tmpDimsNew(i);
            end

            gridOld = cell(1,obj.numberSpatialDims);
            gridNew = cell(1,obj.numberSpatialDims);

            [gridOld{:}] = ndgrid(idxOld{:});
            [gridNew{:}] = ndgrid(idxNew{:});

            stringResize = '';
            stringResize2 = '';
            for j=1:numel(tmpDimsOld)
                stringResize = [stringResize,',gridOld{',num2str(j),'}'];
                stringResize2 = [stringResize2,',gridNew{',num2str(j),'}'];
            end
            stringResize = stringResize(2:end);
            
            %upsample fields
            for j=1:numel(obj.flexboxObj.x)
                obj.flexboxObj.x{j} = reshape(obj.flexboxObj.x{j},tmpDimsOld);
                
                resizedVar = eval(['interpn(',stringResize,',obj.flexboxObj.x{',num2str(j),'}',stringResize2,',''nearest'');']);
                
                check1 = isnan(resizedVar);
                check2 = isinf(resizedVar);
                
                if ((sum(check1(:)) + sum(check2(:))) > 0)
                    error('Isinf + Isnan!!!')
                end
                
                dimNumber = mod(j,numel(tmpDimsOld)) + 1;
                
                obj.flexboxObj.x{j} = resizedVar*tmpDimsNew(dimNumber)/tmpDimsOld(dimNumber);
                obj.flexboxObj.x{j} = obj.flexboxObj.x{j}(:);
                
                obj.flexboxObj.xBar{j} = obj.flexboxObj.x{j};

                %set dims
                obj.flexboxObj.dims{j} = tmpDimsNew;
            end

            for j=1:numel(obj.flexboxObj.y)
                obj.flexboxObj.y{j} = reshape(obj.flexboxObj.y{j},tmpDimsOld);
                
                resizedVar = eval(['interpn(',stringResize,',obj.flexboxObj.y{',num2str(j),'}',stringResize2,',''nearest'');']);
                
                obj.flexboxObj.y{j} = resizedVar(:);
                
                check1 = isnan(resizedVar);
                check2 = isinf(resizedVar);
                
                if ((sum(check1(:)) + sum(check2(:))) > 0)
                    error('Isinf + Isnan!!!')
                end
            end

            %replace data terms and regularizer
            for j=1:numel(obj.flowTermNumber)
                if (strcmp(obj.dataTerm,'L2'))
                    flowTermObject = L2opticalFlowTerm(1,obj.listImages{level,j},obj.listImages{level,j+1},'discretization',obj.imageDiscretization);
                else %default is L1
                    flowTermObject = L1opticalFlowTerm(1,obj.listImages{level,j},obj.listImages{level,j+1},'discretization',obj.imageDiscretization);
                end

                obj.flexboxObj.duals{obj.flowTermNumber(j)} = flowTermObject;
                
                clear flowTermObject;

                if (obj.doGradientConstancy)
                    for k=1:size(obj.gradConstancyTermNumber,2)
                        if (strcmp(obj.dataTerm,'L2'))
                            gradientConstancyTermObject = L2opticalFlowTerm(1,obj.listImages{level,j},obj.listImages{level,j+1},'termType','gradientConstancy','constancyDimension',k);
                        else %default is L1
                            gradientConstancyTermObject = L1opticalFlowTerm(1,obj.listImages{level,j},obj.listImages{level,j+1},'termType','gradientConstancy','constancyDimension',k);
                        end

                        obj.flexboxObj.duals{obj.gradConstancyTermNumber(j,k)} = gradientConstancyTermObject;
                    end
                    
                    clear gradientConstancyTermObject;
                end
            end

            %adjust regularization
            if (strcmp(obj.regularizerTerm,'L2'))
                regularizerObject = L2gradient(obj.alpha,tmpDimsNew);
            elseif (strcmp(obj.regularizerTerm,'TV'))
                regularizerObject = L1gradientIso(obj.alpha,tmpDimsNew);
            else %default is Huber
                regularizerObject = huberGradient(obj.alpha,tmpDimsNew,0.01);
            end

            for j=1:size(obj.regularizerTermNumber,1)
                for k=1:size(obj.regularizerTermNumber,2)
                    obj.flexboxObj.duals{obj.regularizerTermNumber(j,k)} = regularizerObject;
                end
            end
            
            clear regularizerObject;
            
            for warps = 1:obj.numberOfWarps
                if (obj.doWarping)
                    v = obj.getResult;
                    %warp data terms
                    for j=1:obj.numberImages-1
                        idx = repmat({':'},1,obj.numberSpatialDims + 2);
                        idx{end-1} = j;
                        
                        obj.flexboxObj.duals{obj.flowTermNumber(j)}.warpDataterm(squeeze(v(idx{:})));

                        if (obj.doGradientConstancy)
                            for k=1:obj.numberSpatialDims
                                obj.flexboxObj.duals{obj.gradConstancyTermNumber(j,k)}.warpDataterm(squeeze(v(idx{:})));
                            end
                        end
                    end
                end

                obj.flexboxObj.runAlgorithm;
                
                if (obj.medianFiltering)
                    for j=1:obj.numberImages-1
                        for k=1:obj.numberSpatialDims
                            varNumber = obj.numberSpatialDims*(j-1)+k;
                            if (obj.numberSpatialDims == 2)
                                obj.flexboxObj.x{varNumber} = medfilt2(obj.flexboxObj.getPrimal(varNumber), [5 5],'symmetric');
                            elseif (obj.numberSpatialDims == 3)
                                obj.flexboxObj.x{varNumber} = medfilt3(obj.flexboxObj.getPrimal(varNumber), [5 5 5],'symmetric');
                            else
                                error('Median filtering only exists for dimension 2 and 3. Please turn this option off!');
                            end
                            obj.flexboxObj.x{varNumber} = obj.flexboxObj.x{varNumber}(:);
                        end

                        if (obj.verbose > 1)%on massive verbose
                            clear flowField;
                            
                            if (obj.numberSpatialDims == 2)
                                flowField(:,:,2) = obj.flexboxObj.getPrimal(2*j-1);
                                flowField(:,:,1) = obj.flexboxObj.getPrimal(2*j);

                                figure(300+j);imagesc(flowToColorV2(flowField));axis image;title(['Level #',num2str(i),' : Field #',num2str(j)]);drawnow;
                            elseif (obj.numberSpatialDims == 3)
                                %one slice
                                flowField1 = obj.flexboxObj.getPrimal(3*j-1);
                                flowField2 = obj.flexboxObj.getPrimal(3*j);
                                
                                nSlices = size(flowField1,3);
                                flowField1 = flowField1(:,:,floor(nSlices/2));
                                flowField2 = flowField2(:,:,floor(nSlices/2));
                                
                                figure(300+j);imagesc(flowToColorV2(cat(3,flowField1,flowField2)));axis image;
                                title(['Level #',num2str(i),' : Field #',num2str(j)]);drawnow;
                                
                            end
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
                obj.runLevel(i);
            end
            
            if (obj.verbose > 0)
                disp('Finished');
            end
        end
        
    end
    
end


