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
        temp = zeros(sum(segs{k}(:)),modeNumber);
        for m=1:modeNumber
             temp0 = detectorWave(:,:,m);
             temp(:,m) = temp0(segs{k}==1);
        end
        temp2 = sqrt(sqrtInt(segs{k}==1).^2./(sum(sum(abs(temp).^2))+eps)).*temp;
        for m=1:modeNumber
             temp0 = correctedWave(:,:,m);
             temp0(segs{k}==1) = temp2(:,m);
             correctedWave(:,:,m) =temp0;
        end
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
        temp = zeros(sum(segs{k}(:)),modeNumber);
        for m=1:modeNumber
             temp0 = detectorWave(:,:,m);
             temp(:,m) = temp0(segs{k}==1);
        end
        temp2 = sqrt(sqrtInt(segs{k}==1).^2./(sum(sum(abs(temp).^2))+eps)).*temp;
        for m=1:modeNumber
             temp0 = correctedWave(:,:,m);
             temp0(segs{k}==1) = temp2(:,m);
             correctedWave(:,:,m) =temp0;
        end
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
