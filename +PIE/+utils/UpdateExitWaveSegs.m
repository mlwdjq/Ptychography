function [exitWaveNew,detectorWave,segInt] = UpdateExitWaveSegs(exitWave,sqrtInt,segs,propagator,H,Hm,preShift,dx,wavelength,z)
if nargin <8
    detectorWave = PIE.utils.postPropagate (exitWave,propagator,H,preShift);
    correctedWave = detectorWave;
    for k=1:length(segs)
        correctedWave(segs{k}==1) = sqrt(sqrtInt(segs{k}==1).^2./(sum(abs(detectorWave(segs{k}==1)).^2)+eps)).*detectorWave(segs{k}==1);
    end
    exitWaveNew = PIE.utils.postPropagate (correctedWave,propagator,Hm,preShift);
else
    detectorWave = PIE.utils.Propagate(exitWave,propagator,dx,wavelength,z);
    correctedWave = detectorWave;
    for k=1:length(segs)
        correctedWave(segs{k}==1) = sqrt(sqrtInt(segs{k}==1).^2./(sum(abs(detectorWave(segs{k}==1)).^2)+eps)).*detectorWave(segs{k}==1);
    end
    exitWaveNew = PIE.utils.Propagate(correctedWave,propagator,dx,wavelength,-z);
end

% segments
detectorWave2 = abs(detectorWave).^2;
segInt=zeros(length(segs),1);
for k=1:length(segs)
    seg = logical(segs{k});
    segInt(k) = sum(sum(detectorWave2(seg)));
end
