% reference: Odstr?il, M., Menzel, A., & Guizar-Sicairos, M. (2018).
% Iterative least-squares solver for generalized maximum-likelihood ptychography. Optics express, 26(3), 3108-3123.

function [dObjectRecon,dProbeRecon,tempError,aj,alpha,beta] = sML(dObjectRecon,dProbeRecon,sqrtInt,Iseg,Rpix,N,propagator,...
    H,Hm, gamma, tempError, likelihoodType,regularization,aj,segs)

reconBox = dObjectRecon(Rpix(1)+[1:N],Rpix(2)+[1:N]);
exitWave = reconBox.*dProbeRecon;
detectorWave = PIE.utils.postPropagate (exitWave,propagator,H,1);
% segments
detectorWave2 = abs(detectorWave).^2;
segInt=zeros(length(segs),1);
for k=1:length(segs)
    seg = logical(segs{k});
    segInt(k) = sum(sum(detectorWave2(seg)));
end

% optmization in reciprocal space
switch likelihoodType
    case 'Poisson'
        cosi = 1-Iseg./(segInt+eps);
        gLike = sum(cosi)*detectorWave;
        aj = 1/N^2*sum(segInt-cosi.*(Iseg./(1-aj*cosi)))...
            /sum(cosi.^2.*segInt);
    case 'amplitude'
        sigma =0.5; % only poisson noise
        gLike = 2/(2*sigma)^2*sum(1-sqrt(Iseg)./(sqrt(segInt)+eps)).*detectorWave;
        aj = 1/(2*N^2);
end
detectorWaveOpt = detectorWave - aj*gLike;

% optmization in real space
correctedWave = detectorWaveOpt;
for k=1:length(segs)
    correctedWave(segs{k}==1) = sqrt(sqrtInt(segs{k}==1).^2./sum(abs(detectorWaveOpt(segs{k}==1)).^2)+eps).*detectorWaveOpt(segs{k}==1);
end


exitWaveNew = PIE.utils.postPropagate (correctedWave,propagator,Hm,1);
exitWaveDiff = exitWaveNew - exitWave;
tempProbe = dProbeRecon;
denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
dO = conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
dP = conj(reconBox).*(exitWaveNew-exitWave)./denomP;
C =[sum(real(exitWaveDiff(:).*conj(dO(:).*tempProbe(:))));sum(real(exitWaveDiff(:).*conj(dP(:).*reconBox(:))))];
s = dO(:).*tempProbe(:).*conj(dP(:).*reconBox(:))*0;
A = [sum(abs(dO(:).*tempProbe(:)).^2)+regularization ,sum(s) ; ...
    sum(conj(s)),sum(abs(dP(:).*reconBox(:)).^2)+regularization ];
x = pinv(A)*C;
alpha = x(1);
beta = x(2);

newReconBox = reconBox + alpha*dO;
dProbeRecon = dProbeRecon + beta*dP;
dObjectRecon(Rpix(1)+[1:N],Rpix(2)+[1:N]) = newReconBox;
tempError = tempError + abs(sqrt(Iseg)-sqrt(segInt)).^2;
