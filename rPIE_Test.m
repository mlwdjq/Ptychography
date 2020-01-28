%% this script is used for test the rPIE algorithm
% reference: Rodenburg, J., & Maiden, A. (2019). Ptychography.
% In Springer Handbook of Microscopy (pp. 2-2). Springer, Cham.
addpath(genpath('D:\OneDrive\Ptychography\code\mip-master'));
addpath('E:\matlab\WZhu_toolbox');

%% parameters setting
lambda_um =13.5e-3; % wavelength
N = 156; % grid number of calculation box in object plane
dc_um = 8; % detector pixel pitch
propagator = 'angular spectrum';
alpha = 1; % weight factor for updating object
beta = 1; % weight factor for updating probe
z_mm = 200;
z_um = z_mm*1000; % distance from object plane to detector plane
iteration = 1000; % maximum iteration times
samplingFactor_det = lambda_um.*z_um/(dc_um*N*dc_um);
if samplingFactor_det>1
    fprintf('Please adjust configurations for propagation sampling\n');
    return;
end
if strcmp(propagator,'fourier')
    do_um = lambda_um*z_um/N/dc_um;
else
    do_um = dc_um; % object pixel pitch
end
scanSteps = 4; % scanning steps scanSteps*scanSteps
Rmin_um = 0; % minimum shift
Rmax_um = 400; % maximum shift
[Rx_um,Ry_um] = meshgrid(linspace(Rmin_um,Rmax_um,scanSteps));
Rx_um(2:2:end,:) =fliplr(Rx_um(2:2:end,:));
Rx_um=Rx_um';
Ry_um=Ry_um';
Rpix = round(([Ry_um(:),Rx_um(:)]-Rmin_um)/do_um);
K = max(Rpix(:,1))-min(Rpix(:,1))+N;
L = max(Rpix(:,2))-min(Rpix(:,2))+N; % size of object [K,L]
xp_um = [1:N]*do_um; % object coordinates
xo_um = [1:L]*do_um; % object coordinates
yo_um = [1:K]*do_um; % object coordinates
xc_um = [1:N]*dc_um; % detector coordinates
xp_mm = xp_um/1000;
xo_mm = xo_um/1000;
yo_mm = yo_um/1000;
xc_mm = xc_um/1000;

%% initial probe
probeType = 'defocus';
switch probeType
    case 'defocus'
        df_um = -50e3; % negative sign corresponds convergent
        samplingFactor_obj = lambda_um.*z_um/(N*dc_um*do_um);
        Rc_um = (N-20)*dc_um/2;
        Rprobe_um = abs(df_um)/(z_um+df_um)*Rc_um;
        [n1,n2]=meshgrid(1:N);
        n1 = n1-N/2-1;
        n2 = n2-N/2-1;
        probe = ifftshift(ifft2(ifftshift(pinhole(round(Rc_um/dc_um),N,N).*...
            exp(-1i*pi*df_um*dc_um^2/lambda_um/z_um^2*(n1.^2+n2.^2)))));
    case 'plane'
        Rprobe_um = do_um*N*0.15;
        dp_um = 50000; % propagation distance
        samplingFactor_obj = lambda_um.*dp_um/(N*do_um*do_um);
        probe = Propagate (pinhole(round(2*Rprobe_um/do_um),N,N),'angular spectrum',...
            do_um,lambda_um,dp_um);
end
if samplingFactor_obj>1
    fprintf('Please adjust configurations for propagation sampling\n');
    return;
end
overlap = overlapRatio(Rprobe_um,Rmax_um/scanSteps); % overlap ratio of two circles
if overlap<0.6 % check overlap
    fprintf('Overlap ratio(%0.1f%%) is less than 60%%. Please adjust configurations\n',overlap*100);
    return;
end
% figure(2),imagesc(xp_um,xp_um,atan2(imag(probe),real(probe)));
% xlabel('mm'),ylabel('mm');axis tight equal
% figure(3),imagesc(xp_mm,xp_mm,abs(probe)),xlabel('mm'),ylabel('mm');axis tight equal

%% initial object
object_amp = mat2gray(imread('cameraman.tif'))*0.8+0.2; % object amplitude
object_amp = object_amp(1:K,1:L);
object_pha=mat2gray(imread('pears.png'));
object_pha= object_pha(1:K,1:L,1);
object_pha = (object_pha-0.5)*2*pi; % object phase
object = object_amp.*exp(1i*object_pha);
% figure(2),imagesc(xo_mm,yo_mm,abs(object_amp)),xlabel('mm'),ylabel('mm');axis tight equal

%% simulate diffracted patterns
Em = zeros(N,N,scanSteps^2);
for j = 1:scanSteps^2
    reconBox = object(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
    exitWave = reconBox.*probe;
    Em(:,:,j) = Propagate (exitWave ,propagator,do_um,lambda_um,z_um);
    measurements = abs(Em).^2;
%     figure(2),imagesc(xc_mm,xc_mm,measurements(:,:,j)),xlabel('mm'),ylabel('mm');
%     axis tight equal;drawnow;
end
Is = sum(measurements,3); % total intensity on detector
fprintf('diffraction simulation finished\n');

object_sim = object;
probe_sim = probe;

%% define initial guesses for probe and object
object = ones(K,L);
% probe = ones(N);
% probe = pinhole(round(Rprobe_um/do_um),N,N);

%% reconstruction
method = 'rPIE';
errors = zeros(iteration,1);
gamma = 0.2; % weight factor for rPIE, 1 for ePIE
delta = 0.1; % weight factor for RAAR, 1 for DM
for i = 1:iteration % doing iteration
    tempError = 0;
    switch method
        case 'rPIE'
            for j =1:scanSteps^2 % scanning solution
                reconBox = object(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
                % figure(2),imagesc(xo_mm,yo_mm,abs(object)),axis tight equal ;
                % title('object amplitude');drawnow;
                exitWave = reconBox.*probe;
                [exitWaveNew,detectorWave] = UpdateExitWave(exitWave,measurements(:,:,j),...
                    propagator,do_um,lambda_um,z_um);
                tempProbe = probe;
                denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
                newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
                denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
                probe = probe + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
                object(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]) = newReconBox;
                tempError = tempError + abs(sqrt(measurements(:,:,j))-abs(detectorWave)).^2;
            end
        case 'RAAR'
            if i==1 % initial exitWaves
                exitWaves = zeros(N,N,scanSteps^2);
                for j =1:scanSteps^2
                    reconBox = object(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
                    exitWaves(:,:,j) = probe.*reconBox;
                end
            end
            parfor j =1:scanSteps^2 % batch scanning solution
                reconBox = object(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
                waveToPropagate = 2*probe.*reconBox-exitWaves(:,:,j);
                [exitWaveNew,detectorWave] = UpdateExitWave(waveToPropagate,measurements(:,:,j),...
                    propagator,do_um,lambda_um,z_um);
                tempError = tempError + abs(sqrt(measurements(:,:,j))-abs(detectorWave)).^2;
                exitWaves(:,:,j) = delta*(exitWaves(:,:,j) + exitWaveNew) +(1-2*delta)*probe.*reconBox;
            end
            probe = BatchProbeUpdate(exitWaves,probe,object,Rpix,beta);
            object = BatchObjectUpdate(exitWaves,probe,object,Rpix,alpha);
    end
    
    % error evaluation
    errors(i) = sum(tempError(:))/sum(Is(:));
    
    % intermedium object and probe
    Eobj_amp = abs(object);
    Eobj_pha = atan2(imag(object),real(object));
    Epro_amp = abs(probe);
    Epro_pha = atan2(imag(probe),real(probe));
    figure(3),subplot(221),imagesc(xo_mm,yo_mm,Eobj_amp),axis tight equal ;title('Object amplitude');
    xlabel('mm'),  ylabel('mm')
    subplot(222),imagesc(xo_mm,yo_mm,Eobj_pha),axis tight equal ;title('Object phase');xlabel('mm'),  ylabel('mm')
    subplot(223),imagesc(xp_mm,xp_mm,Epro_amp),axis tight equal ;title('Probe amplitude');xlabel('mm'),  ylabel('mm')
    subplot(224),imagesc(xp_mm,xp_mm,Epro_pha),axis tight equal ;title('probe phase');xlabel('mm'),  ylabel('mm')
    drawnow;
    fprintf('%d iterations finished,residual error: %0.5f\n',i,errors(i));
    if i>1&&((abs(errors(i)-errors(i-1))<1e-7)||errors(i)<1e-4)
        break;
    end
end

