% Removes spherical term from wavefront by minimizing rms residual


function [dWaveOut, dFVal] = removeSphere(dWaveIn, NA)

[sr, sc] = size(dWaveIn);
if (sr ~= sc)
    error('Input wave must be square');
end

%dMsk = pinhole(sr);
dMsk=pinhole(sr);
idx = linspace(-1, 1, sr);
[dX, dY] = meshgrid(idx);
dZ = 1/(tan(asin(NA)));

hSphere =  sqrt(dX.^2 + dY.^2 + dZ.^2);
dIdx = dMsk > 0;
hObjective = @(dC) std((dWaveIn(dIdx) - dC*hSphere(dIdx).*dMsk(dIdx)));

[dCOpt, dFVal] = fminsearch(hObjective, 0);

dWaveOut = dWaveIn - dCOpt*hSphere;