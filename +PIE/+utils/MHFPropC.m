% Ryan Miyakawa 8/2013
%
% This function requires the compiled MATLAB executable MHFPropCmex to run.
% This file is mex-compiled from the source: MHFPropCmex.cpp
%
% function output =
%              MHFPropC(aperture, inCoordsX, inCoordsY, inCoordsZ,
%              outCoordX, outCoordY, outCoordZ, wavelength_nm)
%               
%
% Full propagation to customized coordinates from customized coordinates
% using MEX for speed
%
%



function output = MHFPropC(aperture, inCoordsX, inCoordsY, inCoordsZ, outCoordX, outCoordY, outCoordZ, wavelength)

if nargin ~= 8
    fprintf('**Bad input**\n\nEnter as:\nMHFPropC(ap, inX, inY, inZ, outX, outY, outZ, lam_nm)  or\n');
    fprintf('MHFPropC(ap, inX, inY, [], outX, outY, z_um, lam_nm)  or\n');
    fprintf('MHFPropC(ap, inX, [], [], outX, [], z_um, lam_nm)  or\n');
    error('Bad Input');
end



% Allow Y coords to be []
if isempty(inCoordsY)
    inCoordsY = zeros(size(inCoordsX));
end
if isempty(outCoordY)
    outCoordY = zeros(size(outCoordX));
end

if (isempty(inCoordsZ) && length(outCoordZ(:)) == 1)
    inCoordsZ = zeros(size(aperture));
    outCoordZ = outCoordZ * ones(size(outCoordX));
end

% Check sampling:
HFSampling(aperture, inCoordsX, inCoordsY, inCoordsZ, outCoordX, outCoordY, outCoordZ, wavelength * 1000);

% Only pass nonzero components:
idx = find(aperture);
tic;
output = MHFPropCmex(aperture(idx), inCoordsX(idx), inCoordsY(idx), inCoordsZ(idx), ...
                        outCoordX, outCoordY, outCoordZ, wavelength * 1000);
fprintf('Propagation took %s\n', s2f(toc));