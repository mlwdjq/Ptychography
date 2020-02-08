function [dObjectRecon,dProbeRecon,tempError,Cpix,rot_rad,scale] = pcPIE(dObjectRecon,...
    dProbeRecon,sqrtInt,tempError,alpha,beta,gamma,propagator,H,Hm,N,KL,Rpix,...
    Cpix,rot_rad,scale, nC,rCpix0,maxRot_rad0,maxScale0,decreaseFactor,centralPix,correctMethod)
detectorWaveTry = zeros(N,N,nC);
posError = zeros(nC,1);
switch correctMethod
    case 'pcPIE' % annealing algorithm
        rCpix = rCpix0*decreaseFactor;
        randPix= [[0,0];round(rCpix*(rand(nC-1,2)*2-1))]; % generate rand searching pixel within the radius rCpix
        RpixTry = Rpix + Cpix + randPix;
        RpixTry( RpixTry <= 0) = 0;
        RpixTry( RpixTry+N>KL) = 0;
        % find the best position by comparing the diffraction intensity
        for m =1:nC
            reconBox = dObjectRecon(RpixTry(m,1)+[1:N],RpixTry(m,2)+[1:N]);
            exitWave = reconBox.*dProbeRecon;
            detectorWaveTry(:,:,m) = PIE.utils.postPropagate (exitWave,propagator,H,1);
            temp = (abs(detectorWaveTry(:,:,m))-sqrtInt).^2;
            posError(m) = sum(temp(:));
        end
        [~,m0] = min(posError);
        if all(abs(Cpix+randPix(m0,:))<=rCpix0)
            Cpix = Cpix+randPix(m0,:);
        else
            Cpix = [0,0];
        end
    case 'e-pcPIE' % extend pcPIE,which corect additional rotation and scale
        rCpix = rCpix0*decreaseFactor;
        maxRot_rad =maxRot_rad0*decreaseFactor;
        maxScale = maxScale0*decreaseFactor;
        
        randPix= [[0,0];rCpix*(rand(nC-1,2)*2-1)]; % generate rand searching pixel within the radius rCpix
        randRot_rad= [0;maxRot_rad*(rand(nC-1,1)*2-1)]; % generate rand searching rotation
        randScale= [0;maxScale*(rand(nC-1,1)*2-1)]; % generate rand searching rotation
        rotPix(:,2) =  (Rpix(2)-centralPix)*cos(rot_rad + randRot_rad)+...
            (Rpix(1)-centralPix)*sin(rot_rad + randRot_rad)-Rpix(2)+centralPix;
        rotPix(:,1) =  (Rpix(1)-centralPix)*cos(rot_rad + randRot_rad)-...
            (Rpix(2)-centralPix)*sin(rot_rad + randRot_rad)-Rpix(1)+centralPix;
        scalePix = (Rpix-centralPix).*(scale + randScale);
        
        RpixTry = Rpix + round(Cpix+ randPix + rotPix+ scalePix);
        RpixTry( RpixTry <= 0) = 0;
        RpixTry( RpixTry+N>KL) = 0;
        % find the best position by comparing the diffraction intensity
        for m =1:nC
            reconBox = dObjectRecon(RpixTry(m,1)+[1:N],RpixTry(m,2)+[1:N]);
            exitWave = reconBox.*dProbeRecon;
            detectorWaveTry(:,:,m) = PIE.utils.postPropagate (exitWave,propagator,H,1);
            temp = (abs(detectorWaveTry(:,:,m))-sqrtInt).^2;
            posError(m) = sum(temp(:));
        end
        [~,m0] = min(posError);
        
        % restart searching
        if all(abs(Cpix+randPix(m0,:))<=rCpix0)
            Cpix = Cpix+randPix(m0,:);
        else % bad serach
            Cpix = [0,0];
        end
        if abs(rot_rad + randRot_rad(m0))<= maxRot_rad0
            rot_rad = rot_rad + randRot_rad(m0);
        else
            rot_rad =0;
        end
        if abs(scale + randScale(m0))<=maxScale0
            scale = scale + randScale(m0);
        else
            scale =0;
        end
        
end
% update probe and object using new position
detectorWave=detectorWaveTry(:,:,m0);
reconBox = dObjectRecon(RpixTry(m0,1)+[1:N],RpixTry(m0,2)+[1:N]); % get reconstruction region
exitWave = reconBox.*dProbeRecon; % calculate exit wave
correctedWave = sqrtInt.*detectorWave./(abs(detectorWave)+eps); % correct diffraction intensity
exitWaveNew = PIE.utils.postPropagate (correctedWave,propagator,Hm,1); % update exitwave
% update probe and object
tempProbe = dProbeRecon;
denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
dProbeRecon = dProbeRecon + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;% update object
dObjectRecon(RpixTry(m0,1)+[1:N],RpixTry(m0,2)+[1:N]) = newReconBox;% update probe
tempError = tempError + abs(sqrtInt-abs(detectorWave)).^2;