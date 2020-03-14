function [exitWaveNew,detectorWave] = UpdateMultimodeExitWave3D(exitWave,sqrtInt,modeNumber,propagator,H,Hm,preShift,dc_um,lambda_um,Rpix,dx,z)
[i,j,k] = size(exitWave);
detectorWave = zeros(i,j,k);
exitWaveNew = zeros(i,j,k);
if nargin <11
    for m=1:modeNumber
        detectorWave(:,:,m) = PIE.utils.postPropagate (exitWave(:,:,m),propagator,H,preShift);
        detectorWave(:,:,m) = PIE.utils.Propagate(detectorWave(:,:,m),'angular spectrum',dc_um,lambda_um,-Rpix(3)*1000);
    end
    correctedWave = sqrtInt./(sqrt(sum(abs(detectorWave).^2,3))+eps).*detectorWave;
    for m=1:modeNumber
        correctedWave(:,:,m) = PIE.utils.Propagate(correctedWave(:,:,m),'angular spectrum',dc_um,lambda_um,Rpix(3)*1000);
        exitWaveNew(:,:,m) = PIE.utils.postPropagate (correctedWave(:,:,m),propagator,Hm,preShift);
    end
else
    for m=1:modeNumber
        detectorWave(:,:,m)  = PIE.utils.Propagate(exitWave(:,:,m),propagator,dx,lambda_um,z);
        detectorWave(:,:,m) = PIE.utils.Propagate(detectorWave(:,:,m),'angular spectrum',dc_um,lambda_um,-Rpix(3)*1000);
    end
    correctedWave = sqrtInt.*detectorWave./(sqrt(sum(abs(detectorWave).^2,3))+eps);
    for m=1:modeNumber
        correctedWave(:,:,m) = PIE.utils.Propagate(correctedWave(:,:,m),'angular spectrum',dc_um,lambda_um,Rpix(3)*1000);
        exitWaveNew(:,:,m) = PIE.utils.Propagate(correctedWave(:,:,m),propagator,dx,lambda_um,-z);
    end
end