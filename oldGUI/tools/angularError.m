function [ error ] = angularError(  v1,vGT )
    %Computes the angular error between two velocity fields v1 and v2 per
    %pixel
    
    %don't compare for unknown flow
    v1(vGT > 1e5) = 0;
    vGT(vGT > 1e5) = 0;
    
    %project v1 and v2 into higher dimensional space and normalize
    v1(:,:,3) = ones(size(v1,1),size(v1,2));
    vGT(:,:,3) = ones(size(v1,1),size(v1,2));
    
    normV1 = sqrt(v1(:,:,1).^2 + v1(:,:,2).^2 + 1);
    normV2 = sqrt(vGT(:,:,1).^2 + vGT(:,:,2).^2 + 1);

    %v1 = v1 ./ cat(3,normV1,normV1,normV1);
    %v2 = v2 ./ cat(3,normV2,normV2,normV2);
	
	cross = v1(:,:,1).*vGT(:,:,1) + v1(:,:,2).*vGT(:,:,2) + v1(:,:,3).*vGT(:,:,3);

    error = acos(cross ./ (normV1.*normV2));
    error = 2*sum(error(:)) / numel(v1);
end

