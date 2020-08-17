function [exitWaveNew,detectorWave,segInt] = UpdateMultimodeExitWaveSegs(exitWave,sqrtInt,modeNumber,segs,propagator,H,Hm,preShift,dx,wavelength,z)
[i,j,k] = size(exitWave);
detectorWave = zeros(i,j,k);
exitWaveNew = zeros(i,j,k);
if nargin <9
    for m=1:modeNumber
        detectorWave(:,:,m) = PIE.utils.postPropagate (exitWave(:,:,m),propagator,H{m},preShift);
    end
    correctedWave = detectorWave;
    for k=1:length(segs)
        flag = repmat(segs{k}==1,1,1,modeNumber);
        correctedWave(flag) = sqrt(repmat(sqrtInt(segs{k}==1).^2,modeNumber,1)./...
            (sum(abs(detectorWave(flag)).^2)+eps)).*detectorWave(flag);
    end
    for m=1:modeNumber
        exitWaveNew(:,:,m) = PIE.utils.postPropagate (correctedWave(:,:,m),propagator,Hm{m},preShift);
    end
else
    for m=1:modeNumber
        detectorWave = PIE.utils.Propagate(exitWave(:,:,m),propagator,dx,wavelength(m),z);
    end
    correctedWave = detectorWave;
    for k=1:length(segs)
        correctedWave(repmat(segs{k}==1,1,1,modeNumber)) = sqrt(repmat(sqrtInt(segs{k}==1).^2,modeNumber,1)./...
            (sum(abs(detectorWave(repmat(segs{k}==1,1,1,modeNumber))).^2)+eps)).*detectorWave(repmat(segs{k}==1,1,1,modeNumber));
    end
    for m=1:modeNumber
        exitWaveNew(:,:,m) = PIE.utils.Propagate(correctedWave(:,:,m),propagator,dx,wavelength(m),-z);
    end
end

% segments
detectorWave2 = abs(detectorWave).^2;
segInt=zeros(length(segs),1);
for k=1:length(segs)
    seg = logical(segs{k});
    for m = 1:modeNumber
        temp = detectorWave2(:,:,m);
        segInt(k) = segInt(k) + sum(sum(temp(seg)));
    end
end
