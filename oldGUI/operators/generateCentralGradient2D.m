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
function [ Dx2D,Dy2D ] = generateCentralGradient2D( dims,stepsize )
    ex = ones(dims(1), 1);
    ey = ones(dims(2), 1);

    Dx = 1/(2*stepsize(1)) .* spdiags([-ex ex], [-1,1], dims(1), dims(1));
    Dy = 1/(2*stepsize(2)) .* spdiags([-ey ey], [-1,1], dims(2), dims(2)); 
    
    %Neuman Boundary
    Dx(1,:) = 0;
    Dy(1,:) = 0;
    Dx(end,:) = 0;
    Dy(end,:) = 0;
    
    Dx2D = kron(speye(dims(2)),Dx);
    Dy2D = kron(Dy,speye(dims(1)));
end