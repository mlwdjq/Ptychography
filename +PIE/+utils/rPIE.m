function [dObjectRecon,dProbeRecon,tempError] = rPIE(dObjectRecon,dProbeRecon,sqrtInt,Rpix,N,propagator,...
    H,Hm,alpha,beta, gamma, tempError, modeNumber)
if modeNumber<=1
    reconBox = dObjectRecon(Rpix(1)+[1:N],Rpix(2)+[1:N]);
    exitWave = reconBox.*dProbeRecon;
    [exitWaveNew,detectorWave] = PIE.utils.UpdateExitWave(exitWave,sqrtInt,...
        propagator,H,Hm,1);
    tempProbe = dProbeRecon;
    denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
    newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
    denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
    dProbeRecon = dProbeRecon + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
    dObjectRecon(Rpix(1)+[1:N],Rpix(2)+[1:N]) = newReconBox;
    tempError = tempError + abs(sqrtInt-abs(detectorWave)).^2;
else
    reconBox = dObjectRecon(Rpix(1)+[1:N],Rpix(2)+[1:N]);
    exitWave = reconBox.*dProbeRecon;
    [exitWaveNew,detectorWave] = PIE.utils.UpdateMultimodeExitWave(exitWave,sqrtInt,modeNumber,...
        propagator,H,Hm,1);
    tempProbe = dProbeRecon;
    denomO = gamma*max(max(abs(tempProbe).^2)) + (1-gamma)*abs(tempProbe).^2;
    newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
    denomP = gamma*max(max(abs(reconBox).^2)) + (1-gamma).*abs(reconBox).^2;
    dProbeRecon = dProbeRecon + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
    dObjectRecon(Rpix(1)+[1:N],Rpix(2)+[1:N],:) = newReconBox;
    tempError = tempError + abs(sqrtInt-sqrt(sum(abs(detectorWave).^2,3))).^2;
end