function [dObjectRecon,dProbeRecon,tempError] = sPIE(dObjectRecon,dProbeRecon,sqrtInt,Iseg,Rpix,N,propagator,...
    H,Hm,alpha,beta, gamma, tempError,modeNumber,segs,do_um,lambda_um)
if modeNumber<=1
    reconBox = dObjectRecon(Rpix(1,1)+[1:N],Rpix(1,2)+[1:N]);
    if size(Rpix,2)==3 % 3D scanning
        dProbeRecon = PIE.utils.Propagate(dProbeRecon,'angular spectrum',do_um(1),lambda_um(1),Rpix(1,3));
    end
    exitWave = reconBox.*dProbeRecon;
    [exitWaveNew,~,segInt] = PIE.utils.UpdateExitWaveSegs(exitWave,sqrtInt,segs,...
        propagator,H{1},Hm{1},1);
    tempProbe = dProbeRecon;
    denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
    newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
    denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
    dProbeRecon = dProbeRecon + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
    if size(Rpix,2)==3 % 3D scanning
        dProbeRecon = PIE.utils.Propagate(dProbeRecon,'angular spectrum',do_um(1),lambda_um(1),-Rpix(1,3));
    end
    dObjectRecon(Rpix(1,1)+[1:N],Rpix(1,2)+[1:N]) = newReconBox;
    tempError = tempError + abs(sqrt(Iseg)-sqrt(segInt)).^2;
else
    reconBox =dProbeRecon;
    for m = 1:modeNumber
        reconBox(:,:,m) = dObjectRecon(Rpix(m,1)+[1:N],Rpix(m,2)+[1:N],m);
    end
    if size(Rpix,2)==3 % 3D scanning
        for m = 1:modeNumber
            dProbeRecon(:,:,m) = PIE.utils.Propagate(dProbeRecon(:,:,m),'angular spectrum',do_um(m),lambda_um(m),Rpix(m,3));
        end
    end
    exitWave = reconBox.*dProbeRecon;
    [exitWaveNew,~,segInt] = PIE.utils.UpdateMultimodeExitWaveSegs(exitWave,sqrtInt,modeNumber,segs,...
        propagator,H,Hm,1);
    tempProbe = dProbeRecon;
    denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
    newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
    denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
    dProbeRecon = dProbeRecon + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
    if size(Rpix,2)==3 % 3D scanning
        for m = 1:modeNumber
            dProbeRecon(:,:,m) = PIE.utils.Propagate(dProbeRecon(:,:,m),'angular spectrum',do_um(m),lambda_um(m),-Rpix(m,3));
        end
    end
    for m =1:modeNumber
        dObjectRecon(Rpix(m,1)+[1:N],Rpix(m,2)+[1:N],m) = newReconBox(:,:,m);
    end
    tempError = tempError + abs(sqrt(Iseg)-sqrt(segInt)).^2;
end