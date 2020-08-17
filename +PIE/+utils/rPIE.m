function [dObjectRecon,dProbeRecon,tempError] = rPIE(dObjectRecon,dProbeRecon,sqrtInt,Rpix,N,propagator,...
    H,Hm,alpha,beta, gamma, tempError, modeNumber,do_um,lambda_um,FP,CTF)
if FP==0 || nargin<17
    if modeNumber<=1
        reconBox = dObjectRecon(Rpix(1,1)+[1:N],Rpix(1,2)+[1:N]);
        if size(Rpix,2)==3 % 3D scanning
            dProbeRecon = PIE.utils.Propagate(dProbeRecon,'angular spectrum',do_um(1),lambda_um(1),Rpix(1,3));
        end
        exitWave = reconBox.*dProbeRecon;
        [exitWaveNew,detectorWave] = PIE.utils.UpdateExitWave(exitWave,sqrtInt,...
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
        tempError = tempError + abs(sqrtInt-abs(detectorWave)).^2;
    else
        reconBox = dProbeRecon;
        for m=1:modeNumber
            reconBox(:,:,m) = dObjectRecon(Rpix(m,1)+[1:N],Rpix(m,2)+[1:N],m);
        end
        if size(Rpix,2)==3 % 3D scanning
            for m=1:modeNumber
                dProbeRecon(:,:,m) = PIE.utils.Propagate(dProbeRecon(:,:,m),'angular spectrum',do_um(m),lambda_um(m),Rpix(m,3));
            end
        end
        exitWave = reconBox.*dProbeRecon;
        [exitWaveNew,detectorWave] = PIE.utils.UpdateMultimodeExitWave(exitWave,sqrtInt,modeNumber,...
            propagator,H,Hm,1);
        tempProbe = dProbeRecon;
        denomO = gamma*max(max(abs(tempProbe).^2)) + (1-gamma)*abs(tempProbe).^2;
        newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
        denomP = gamma*max(max(abs(reconBox).^2)) + (1-gamma).*abs(reconBox).^2;
        dProbeRecon = dProbeRecon + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
        for m=1:modeNumber
            dObjectRecon(Rpix(m,1)+[1:N],Rpix(m,2)+[1:N],m) = newReconBox(:,:,m);
        end
        if size(Rpix,2)==3 % 3D scanning
            for m=1:modeNumber
                dProbeRecon(:,:,m) = PIE.utils.Propagate(dProbeRecon(:,:,m),'angular spectrum',do_um(m),lambda_um(m),-Rpix(m,3));
            end
        end
        tempError = tempError + abs(sqrtInt-sqrt(sum(abs(detectorWave).^2,3))).^2;
    end
else
    if modeNumber<=1
        reconBox = dObjectRecon(Rpix(1,1)+[1:N],Rpix(1,2)+[1:N]);
        if size(Rpix,2)==3 % 3D scanning
            dProbeRecon = PIE.utils.Propagate(dProbeRecon,'angular spectrum',do_um(1),lambda_um(1),Rpix(1,3));
        end
        exitWave = reconBox.* dProbeRecon.*CTF;
        [exitWaveNew,detectorWave] = PIE.utils.UpdateExitWave(exitWave,sqrtInt,...
            propagator,H{1},Hm{1},1);
         exitWaveNew = exitWaveNew.*CTF;
        tempProbe = dProbeRecon;
        denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
        newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
        denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
        dProbeRecon = dProbeRecon + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
        if size(Rpix,2)==3 % 3D scanning
            dProbeRecon = PIE.utils.Propagate(dProbeRecon,'angular spectrum',do_um(1),lambda_um(1),-Rpix(1,3));
        end
        dObjectRecon(Rpix(1,1)+[1:N],Rpix(1,2)+[1:N]) = newReconBox;
        tempError = tempError + abs(sqrtInt-abs(detectorWave)).^2;
    else
        reconBox =dProbeRecon;
        for m=1:modeNumber
            reconBox(:,:,m) = dObjectRecon(Rpix(m,1)+[1:N],Rpix(m,2)+[1:N],m);
        end
        if size(Rpix,2)==3 % 3D scanning
            for m=1:modeNumber
                dProbeRecon(:,:,m) = PIE.utils.Propagate(dProbeRecon(:,:,m),'angular spectrum',do_um(m),lambda_um(m),Rpix(m,3));
            end
        end
        exitWave = reconBox.*dProbeRecon.*CTF;
        [exitWaveNew,detectorWave] = PIE.utils.UpdateMultimodeExitWave(exitWave,sqrtInt,modeNumber,...
            propagator,H,Hm,1);
        exitWaveNew = exitWaveNew.*CTF;
        tempProbe = dProbeRecon;
        denomO = gamma*max(max(abs(tempProbe).^2)) + (1-gamma)*abs(tempProbe).^2;
        newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
        denomP = gamma*max(max(abs(reconBox).^2)) + (1-gamma).*abs(reconBox).^2;
        dProbeRecon = dProbeRecon + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
        if size(Rpix,2)==3 % 3D scanning
            for m=1:modeNumber
                dProbeRecon(:,:,m) = PIE.utils.Propagate(dProbeRecon(:,:,m),'angular spectrum',do_um(m),lambda_um(m),-Rpix(m,3));
            end
        end
        for m=1:modeNumber
            dObjectRecon(Rpix(m,1)+[1:N],Rpix(m,2)+[1:N],m) = newReconBox(:,:,m);
        end
        tempError = tempError + abs(sqrtInt-sqrt(sum(abs(detectorWave).^2,3))).^2;
    end
end