
% function output =
%               HFPropC(aperture, inCoordsX, inCoordsY, inCoordsZ,
%                       outCoordsX, outCoordsY, outCoordsZ, lambda)
%
%
% Full 2D propagation to customized coordinates from customized coordinates



function output = HFPropC(aperture, inCoordsX, inCoordsY, inCoordsZ, outCoordsX, outCoordsY, outCoordsZ, lambda)

if isempty(inCoordsY)
    inCoordsY = zeros(size(inCoordsX));
end
if isempty(outCoordsY)
    outCoordsY = zeros(size(outCoordsX));
end

if (isempty(inCoordsZ) && length(outCoordsZ(:)) == 1)
    inCoordsZ = zeros(size(aperture));
    outCoordsZ = outCoordsZ * ones(size(outCoordsX));
end

[oSr, oSc]  = size(outCoordsX);
nInElm      = length(aperture(:));

% Vectorize
aperture  = aperture(:);
inCoordsX = inCoordsX(:);
inCoordsY = inCoordsY(:);
inCoordsZ = inCoordsZ(:);
outCoordsX = outCoordsX(:);
outCoordsY = outCoordsY(:);
outCoordsZ = outCoordsZ(:);

h = waitbar(0,'Propagating...   0%  ETA: calculating...');
t0 = clock;

output = zeros(size(outCoordsX));
for c = 1:nInElm
    if mod(c, 50) == 0
        progress(c/nInElm, t0, h, 'Propagating...');
    end
    
    if aperture(c) ~= 0
        R       = sqrt( (outCoordsX - inCoordsX(c)).^2 + (outCoordsY - inCoordsY(c)).^2  + (outCoordsZ - inCoordsZ(c)).^2  );
        relm    = aperture(c)*cos( 2*pi/lambda * R);
        celm    = aperture(c)*sin( 2*pi/lambda * R);
        output  = output + (relm + 1i*celm)./R.^2;
    end
end

%output = output* z_um/(1i * l_um);

output = reshape(output, oSr, oSc);


close(h);


