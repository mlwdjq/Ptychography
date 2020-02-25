% reference: Odstr?il, M., Menzel, A., & Guizar-Sicairos, M. (2018). 
% Iterative least-squares solver for generalized maximum-likelihood ptychography. Optics express, 26(3), 3108-3123.

function [dObjectRecon,dProbeRecon,tempError,aj,alpha,beta] = ML(dObjectRecon,dProbeRecon,sqrtInt,Rpix,N,propagator,...
    H,Hm, gamma, tempError, likelihoodType,regularization,aj)
reconBox = dObjectRecon(Rpix(1)+[1:N],Rpix(2)+[1:N]);
exitWave = reconBox.*dProbeRecon;
detectorWave = PIE.utils.postPropagate (exitWave,propagator,H,1);

% optmization in reciprocal space
switch likelihoodType
    case 'Poisson'
        cosi = 1-sqrtInt(:).^2./(abs(detectorWave(:)).^2+eps);
        gLike = sum(cosi)*detectorWave;
        aj = 1/N^2*sum(abs(detectorWave(:)).^2-cosi.*(sqrtInt(:).^2./(1-aj*cosi)))...
            /sum(cosi.^2.*abs(detectorWave(:)).^2);
    case 'amplitude'        
        sigma =0.5; % only poisson noise
        gLike = 2/(2*sigma)^2*sum(1-sqrtInt(:)./(abs(detectorWave(:))+eps)).*detectorWave;
        aj = 1/(2*N^2);
end
detectorWaveOpt = detectorWave - aj*gLike;

% optmization in real space
correctedWave = sqrtInt.*detectorWaveOpt./(abs(detectorWaveOpt)+eps);
exitWaveNew = PIE.utils.postPropagate (correctedWave,propagator,Hm,1);
exitWaveDiff = exitWaveNew - exitWave;
tempProbe = dProbeRecon;
denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
dO = conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
dP = conj(reconBox).*(exitWaveNew-exitWave)./denomP;
C =[sum(real(exitWaveDiff(:).*conj(dO(:).*tempProbe(:))));sum(real(exitWaveDiff(:).*conj(dP(:).*reconBox(:))))];
s = dO(:).*tempProbe(:).*conj(dP(:).*reconBox(:));
A = [sum(abs(dO(:).*tempProbe(:)).^2)+regularization ,sum(s) ; ...
    sum(conj(s)),sum(abs(dP(:).*reconBox(:)).^2)+regularization ];
x = pinv(A)*C;
alpha = x(1);
beta = x(2);
% alpha = sum(real(exitWaveDiff(:).*conj(dO(:).*tempProbe(:))))./(sum(abs(dO(:).*tempProbe(:)).^2)+regularization)
% beta = sum(real(exitWaveDiff(:).*conj(dP(:).*reconBox(:))))./(sum(abs(dP(:).*reconBox(:)).^2)+regularization)
newReconBox = reconBox + alpha*dO;
dProbeRecon = dProbeRecon + beta*dP;
dObjectRecon(Rpix(1)+[1:N],Rpix(2)+[1:N]) = newReconBox;
tempError = tempError + abs(sqrtInt-abs(detectorWave)).^2;
