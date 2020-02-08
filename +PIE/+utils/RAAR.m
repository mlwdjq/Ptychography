function [dObjectRecon,dProbeRecon,tempError] = RAAR(dObjectRecon,dProbeRecon,...
    sqrtInt,exitWaves,tempError,alpha,beta,delta,Rpix,N,propagator,H,Hm,scanSteps)
parfor j =1:scanSteps^2 % parallel compute scanning region
    reconBox = dObjectRecon(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
    waveToPropagate = 2*dProbeRecon.*reconBox-exitWaves(:,:,j);
    [exitWaveNew,detectorWave] = PIE.utils.UpdateExitWave(waveToPropagate,sqrtInt(:,:,j),...
        propagator,H,Hm,1);
    tempError = tempError + abs(sqrtInt(:,:,j)-abs(detectorWave)).^2;
    exitWaves(:,:,j) = delta*(exitWaves(:,:,j) + exitWaveNew) +(1-2*delta)*dProbeRecon.*reconBox;
end
dProbeRecon = PIE.utils.BatchProbeUpdate(exitWaves,dProbeRecon,dObjectRecon,Rpix,beta);
dObjectRecon = PIE.utils.BatchObjectUpdate(exitWaves,dProbeRecon,dObjectRecon,Rpix,alpha);