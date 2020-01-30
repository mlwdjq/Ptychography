% Removes spherical term from wavefront by minimizing rms residual


function [dWaveOut, dFVal] = removeSphere1D(dWaveIn, NA, dX)
if nargin~=3
    sr = length(dWaveIn);
    dX = linspace(-1, 1, sr)';
end
dWaveIn=dWaveIn(:);
dX=dX(:);
dZ = 1/(tan(asin(NA)));

hSphere =  sqrt(dX.^2 + dZ.^2);
hObjective = @(dC) std((dWaveIn - dC*hSphere));

[dCOpt, dFVal] = fminsearch(hObjective, 0);

dWaveOut = dWaveIn - dCOpt*hSphere;