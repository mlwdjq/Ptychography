function [dObjectRecon,dProbeRecon,tempError] = sPIE(dObjectRecon,dProbeRecon,sqrtInt,Iseg,Rpix,N,propagator,...
    H,Hm,alpha,beta, gamma, tempError,segs,do_um,lambda_um)
reconBox = dObjectRecon(Rpix(1)+[1:N],Rpix(2)+[1:N]);
if size(Rpix,2)==3 % 3D scanning
    dProbeRecon = PIE.utils.Propagate(dProbeRecon,'angular spectrum',do_um,lambda_um,Rpix(3));
end
exitWave = reconBox.*dProbeRecon;
[exitWaveNew,~,segInt] = PIE.utils.UpdateExitWaveSegs(exitWave,sqrtInt,segs,...
    propagator,H,Hm,1);
tempProbe = dProbeRecon;
denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
dProbeRecon = dProbeRecon + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
if size(Rpix,2)==3 % 3D scanning
    dProbeRecon = PIE.utils.Propagate(dProbeRecon,'angular spectrum',do_um,lambda_um,-Rpix(3));
end
dObjectRecon(Rpix(1)+[1:N],Rpix(2)+[1:N]) = newReconBox;
tempError = tempError + abs(sqrt(Iseg)-sqrt(segInt)).^2;