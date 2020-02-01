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
stepByStepScan = 0;
if stepByStepScan == 1
    Rx_um(2:2:end,:) =fliplr(Rx_um(2:2:end,:));
    Rx_um=Rx_um';
    Ry_um=Ry_um';
end
Rpix = round(([Ry_um(:),Rx_um(:)]-Rmin_um)/do_um);
Rxpix =round((Rx_um-Rmin_um)/do_um);
Rypix =round((Ry_um-Rmin_um)/do_um);
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
probeType = 'plane';
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
            exp(-1i*pi*abs(df_um)*dc_um^2/lambda_um/z_um^2*(n1.^2+n2.^2)))));
%         probes = ifftshift(ifft2(ifftshift(pinhole(round(Rc_um/dc_um),N,N).*...
%             exp(-1i*pi*df_um*dc_um^2/lambda_um/(z_um+10000)^2*(n1.^2+n2.^2)))));
    case 'plane'
        Rprobe_um = do_um*N*0.24;
        dp_um = 50000; % propagation distance
        samplingFactor_obj = lambda_um.*dp_um/(N*do_um*do_um);
        probe = PIE.utils.Propagate (pinhole(round(2*Rprobe_um/do_um),N,N),'angular spectrum',...
            do_um,lambda_um,dp_um);
    case 'aperture'
        Rprobe_um = do_um*N*0.15;
        probe = pinhole(round(2*Rprobe_um/do_um),N,N);
end
probe = single(probe);
if samplingFactor_obj>1
    fprintf('Please adjust configurations for propagation sampling\n');
    return;
end
overlap = PIE.utils.overlapRatio(Rprobe_um,(Rmax_um-Rmin_um)/(scanSteps-1)); % overlap ratio of two circles
if overlap<0.6 % check overlap
    fprintf('Overlap ratio(%0.1f%%) is less than 60%%. Please adjust configurations\n',overlap*100);
%     return;
end
% figure(2),imagesc(xp_um,xp_um,atan2(imag(probe),real(probe)));
% xlabel('mm'),ylabel('mm');axis tight equal
% figure(3),imagesc(xp_mm,xp_mm,abs(probe)),xlabel('mm'),ylabel('mm');axis tight equal

%% initial object
% object_amp = mat2gray(imread('cameraman.tif'))*0.8+0.2; % object amplitude
% % object_amp = ones(K,L);
% object_amp = object_amp(1:K,1:L);
% object_pha=mat2gray(imread('pears.png'));
% object_pha= object_pha(1:K,1:L,1);
% object_pha = (object_pha-0.5)*2*pi; % object phase
% object = object_amp.*exp(1i*object_pha);
% % figure(2),imagesc(xo_mm,yo_mm,abs(object_amp)),xlabel('mm'),ylabel('mm');axis tight equal
% object = single(object);
I = single(flipud(imread('cameraman.tif')));
[m,n] = meshgrid(linspace(0,1,L),linspace(0,1,K));
[sr,sc,~]= size(I);
[p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
object_amp = interp2(p,q,I(:,:,1),m,n);
object_amp = mat2gray(object_amp)*0.8+0.2; % object amplitude
I =single(imread('pears.png'));
[sr,sc,~]= size(I);
[p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
object_pha = interp2(p,q,I(:,:,1),m,n);
object_pha=mat2gray(object_pha);
object_pha = (object_pha-0.5)*2*pi; % object phase
object = object_amp.*exp(1i*object_pha);
%% prepropagate (calculate frequency responds function)
preShift = 1;
H = PIE.utils.prePropagate (probe,propagator,do_um,lambda_um,z_um,preShift);
Hm = PIE.utils.prePropagate (probe,propagator,do_um,lambda_um,-z_um,preShift);
%% simulate diffracted patterns
Em = zeros(N,N,scanSteps^2);
for j = 1:scanSteps^2
    reconBox = object(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
    exitWave = reconBox.*probe;
    %     Em(:,:,j) = Propagate (exitWave ,propagator,do_um,lambda_um,z_um);
    Em(:,:,j) = PIE.utils.postPropagate (exitWave,propagator,H,preShift);
    sqrtInt = single(abs(Em));
    measurements = sqrtInt.^2;
        figure(2),imagesc(xc_mm,xc_mm,sqrtInt(:,:,j)),xlabel('mm'),ylabel('mm');
        axis tight equal;drawnow;
end
Is = sum(measurements,3); % total intensity on detector
fprintf('diffraction simulation finished\n');

object_sim = object;
probe_sim = probe;

%% define initial guesses for probe and object
object = single(ones(K,L));
% probe = ones(N);
% probe = pinhole(round(Rprobe_um/do_um),N,N);

%% reconstruction
method = 'rPIE';
errors = zeros(iteration,1);
alpha = 1; % weight factor for updating object
beta = 0.1; % weight factor for updating probe
gamma = 0.2; % weight factor for rPIE, 1 for ePIE
delta = 0.1; % weight factor for RAAR, 1 for DM
if strcmp(method,'WDD')
    I_set_u_R = reshape(measurements,N,N,scanSteps,scanSteps);
    I_set_R_u= permute(I_set_u_R,[3, 4, 1, 2]); 
    G_set_u_U = single(zeros(N,N,scanSteps,scanSteps));
    G_set_U_u = single(zeros(scanSteps,scanSteps,N,N));
    H_set_r_U = G_set_u_U;
    H_set_U_r = G_set_U_u;
    L_set_r_R = G_set_u_U;
    L_set_R_r = G_set_U_u;
    Xa = G_set_u_U;
    Q = G_set_u_U;
end
for i = 1:iteration % doing iteration
    tempError = 0;
    switch method
        case 'rPIE' % scanning solution
            for j =1:scanSteps^2
                reconBox = object(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
                % figure(2),imagesc(xo_mm,yo_mm,abs(object)),axis tight equal ;
                % title('object amplitude');drawnow;
                exitWave = reconBox.*probe;
                [exitWaveNew,detectorWave] = PIE.utils.UpdateExitWave(exitWave,sqrtInt(:,:,j),...
                    propagator,H,Hm,preShift);
                tempProbe = probe;
                denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
                newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
                denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
                probe = probe + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
                object(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]) = newReconBox;
                tempError = tempError + abs(sqrtInt(:,:,j)-abs(detectorWave)).^2;
            end
        case 'RAAR' % batch scanning solution
            if i==1 % initial exitWaves
                exitWaves = zeros(N,N,scanSteps^2);
                for j =1:scanSteps^2
                    reconBox = object(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
                    exitWaves(:,:,j) = probe.*reconBox;
                end
            end
            parfor j =1:scanSteps^2
                reconBox = object(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
                waveToPropagate = 2*probe.*reconBox-exitWaves(:,:,j);
                [exitWaveNew,detectorWave] = PIE.utils.UpdateExitWave(waveToPropagate,sqrtInt(:,:,j),...
                    propagator,H,Hm,preShift);
                tempError = tempError + abs(sqrtInt(:,:,j)-abs(detectorWave)).^2;
                exitWaves(:,:,j) = delta*(exitWaves(:,:,j) + exitWaveNew) +(1-2*delta)*probe.*reconBox;
            end
            probe = PIE.utils.BatchProbeUpdate(exitWaves,probe,object,Rpix,beta);
            object = PIE.utils.BatchObjectUpdate(exitWaves,probe,object,Rpix,alpha);
        case 'WDD' % Wigner distribution deconvolution 
            % I_set(u,R), G_set(u,U), H_set(r,U), L_set(r,R)
            for m = 1:N
                for n = 1:N
                    G_set_U_u(:,:,m,n) = fftshift(fft2(fftshift(I_set_R_u(:,:,m,n)))); % from R to U
                end
            end
            G_set_u_U = ipermute(G_set_U_u,[3 4 1 2]);
            for m = 1:scanSteps
                for n = 1:scanSteps
                    L_set_r_R(:,:,m,n) = ifftshift(ifft2(ifftshift(I_set_u_R(:,:,m,n)))); % from u to r
                    H_set_r_U(:,:,m,n) = ifftshift(ifft2(ifftshift(G_set_u_U(:,:,m,n)))); % from u to r
                end
            end
            L_set_R_r = permute(L_set_r_R,[3 4 1 2]);
            H_set_U_r = permute(H_set_r_U,[3 4 1 2]);
            for m = 1:scanSteps
                for n = 1:scanSteps
                    temp1 = fftshift(fft2(fftshift(object(Rypix(m,n)+[1:N],Rxpix(m,n)+[1:N])))); % from R to U
                    temp0 = fftshift(fft2(fftshift(object(Rypix(round((scanSteps+1)/2),...
                        round((scanSteps+1)/2))+[1:N],Rxpix(round((scanSteps+1)/2),round((scanSteps+1)/2))+[1:N]))));
                    Xa(:,:,m,n) = ifftshift(ifft2(ifftshift(conj(temp1).*temp0)));%% may change the conj
                end
            end
            Xq = conj(Xa).*H_set_r_U./(abs(Xa).^2+eps);
            for m = 1:scanSteps
                for n = 1:scanSteps
                    Q(:,:,m,n) = fftshift(fft2(fftshift(Xq(:,:,m,n))));%% may change the conj
                end
            end
            for m = 1:scanSteps
                for n = 1:scanSteps
                    Qs(m,n) =Q(N/scanSteps*(m-1)+1,N/scanSteps*(n-1)+1,m,n);%% may change the conj
                end
            end
            q = ifftshift(ifft2(ifftshift(Qs)));
            q_amp =abs(q);
            q_pha =atan2(imag(q),real(q));
            imagesc(q_amp);colorbar;
%             subplot(122),imagesc(q_pha);
%             Probe = fftshift(fft2(fftshift(probe)));
%             temp = permute(G_set_u_U,[1,2,3,4]);
%             imagesc(abs(temp(:,:,10,10)))
            break;
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
    subplot(224),imagesc(xp_mm,xp_mm,Epro_pha),axis tight equal ;title('Probe phase');xlabel('mm'),  ylabel('mm')
    drawnow;
    fprintf('%d iterations finished,residual error: %0.5f\n',i,errors(i));
    if i>1&&((abs(errors(i)-errors(i-1))<1e-7)||errors(i)<1e-4)
        break;
    end
end

