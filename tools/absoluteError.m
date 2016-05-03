function [ error ] = absoluteError( v1,vGT )
    %Computes the absolute error between two velocity fields v1 and v2 per
    %pixel, also known as EndpointError 

    %don't compare for unknown flow
    v1(abs(vGT) > 1e2) = 0;
    vGT(abs(vGT) > 1e2) = 0;
    
    error = sqrt( (v1(:,:,1)-vGT(:,:,1)).^2 + (v1(:,:,2)-vGT(:,:,2)).^2);
    error = error(2:end-1,2:end-1);
    
    figure(701);clf;imagesc(error);colorbar;
    
    
    error = 2*sum(error(:)) / numel(v1);

end

