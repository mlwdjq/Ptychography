function [exitWaveNew,detectorWave] = UpdateMultimodeExitWave(exitWave,sqrtInt,modeNumber,propagator,H,Hm,preShift,dx,wavelength,z)
[i,j,k] = size(exitWave);
detectorWave = zeros(i,j,k);
exitWaveNew = zeros(i,j,k);
if nargin <8
    for m=1:modeNumber
        detectorWave(:,:,m) = PIE.utils.postPropagate (exitWave(:,:,m),propagator,H,preShift);
    end
    correctedWave = sqrtInt.*detectorWave./(sqrt(sum(abs(detectorWave).^2,3))+eps);
    for m=1:modeNumber
        exitWaveNew(:,:,m) = PIE.utils.postPropagate (correctedWave(:,:,m),propagator,Hm,preShift);
    end
else
    for m=1:modeNumber
        detectorWave(:,:,m)  = PIE.utils.Propagate(exitWave(:,:,m),propagator,dx,wavelength,z);
    end
    correctedWave = sqrtInt.*detectorWave./(sqrt(sum(abs(detectorWave).^2,3))+eps);
    for m=1:modeNumber
        exitWaveNew(:,:,m) = PIE.utils.Propagate(correctedWave(:,:,m),propagator,dx,wavelength,-z);
    end
end