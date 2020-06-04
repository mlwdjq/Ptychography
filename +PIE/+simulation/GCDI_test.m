%% this script is used to verify the validity of the GCDI

%% configurations
N = pie.uieRes.get();
detSize_um = pie.uieDetSize.get()*1000;
z_um       = pie.uiez2.get()*1000;
NA         = pie.uieNA.get();
lambda_um   = pie.uieLambda.get()/1000;
propagator = 'fourier';

%% load object
I = single(flipud(imread('cameraman.tif')));
[m,n] = meshgrid(linspace(0,1,N));
[sr,sc,~]= size(I);
[p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
object_amp = interp2(p,q,I(:,:,1),m,n,'nearest');
object_amp = mat2gray(object_amp)*0.8+0.2; % object amplitude
I =single(imread('pears.png'));
[sr,sc,~]= size(I);
[p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
object_pha = interp2(p,q,I(:,:,1),m,n,'nearest');
object_pha=mat2gray(object_pha);
object_pha = (object_pha-0.5)*1*pi; % object phase
object = object_amp.*exp(1i*object_pha);
figure(2),imagesc(object_amp),colorbar;
figure(3),imagesc(object_pha),colorbar;

%% define probe
Rc_um = z_um*tan(asin(NA));

% generate phase modulation
zernCouples = [40:5:80;ones(1,9)]';
Nz = size(zernCouples,1);
%  Nz = 25;
probe = zeros(N,N,Nz);
for j = 1:Nz
    zernCouple = zernCouples(j,:);
    % generate zernike function form based on zernike couples
    zfn         =   PIE.utils.generateZernikeFunction(zernCouple,N,1);
    [x_um,y_um] = meshgrid(linspace(-detSize_um/2,detSize_um/2,N));
    r = sin(atan(sqrt(x_um.^2+y_um.^2)/z_um))/NA;
    th = atan2(y_um,x_um);
     SLM_pha = 2*pi*zfn(r,th);
%     SLM_pha = 2*pi*randn(N);
%     SLM_pha=2*pi*PIE.utils.generateMSFN(0.6,256,0,5)+SLM_pha*0;
    SLM_amp = ones(N);
    s = cos(2*pi*50*x_um/detSize_um)+cos(2*pi*50*y_um/detSize_um);
%     SLM_amp(s<0)=0;
%     figure(3),imagesc(s)
    SLM(:,:,j) = SLM_amp.*exp(1i*SLM_pha);
    probe(:,:,j) = PIE.utils.Propagate (pinhole(round(2*Rc_um/2/pie.dc_um),N,N).*SLM(:,:,j),'fourier',...
        pie.do_um,lambda_um,1);
    probe_amp = abs(probe(:,:,j));
    probe_pha = atan2(imag(probe(:,:,j)),real(probe(:,:,j)));
       figure(2),imagesc(probe_amp),colorbar;pause(0.1)
% figure(3),imagesc(probe_pha),colorbar;
end





%% simulate patterns
sqrtInt = zeros(N,N,Nz);
for j = 1:Nz
    exitWave = object.*probe(:,:,j);
    Em = PIE.utils.postPropagate (exitWave,'fourier',1,1);
    sqrtInt(:,:,j) = single(abs(Em));
    measurements = sqrtInt.^2;
    figure(3),imagesc(sqrtInt(:,:,j)),colorbar;pause(0.1);
end
totalI = sum(measurements(:));


%% reconstruction
dObjectRecon = ones(N);
alpha =0.5;
beta = 0;
gamma =0.2;

for iteration =1:1000
    tempError=0;
    for j = 1:Nz
        reconBox = dObjectRecon;
        dProbeRecon =probe(:,:,j);
        exitWave = reconBox.*dProbeRecon;
        [exitWaveNew,detectorWave] = PIE.utils.UpdateExitWave(exitWave,sqrtInt(:,:,j),...
            propagator,1,-1,1);
        tempProbe = dProbeRecon;
        denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
        newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
        denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
        dProbeRecon = dProbeRecon + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
        dProbeReconPupil = PIE.utils.Propagate (dProbeRecon,'fourier',...
            pie.do_um,lambda_um,-1);
        if j~=Nz
            js=j+1;
        else
            js =1;
        end
        dProbeReconPupilNext = dProbeReconPupil.*conj(SLM(:,:,j))./(abs(SLM(:,:,j)).^2+eps).*SLM(:,:,js);
        probe(:,:,js) = PIE.utils.Propagate (dProbeReconPupilNext,'fourier',...
            pie.do_um,lambda_um,1);
        dObjectRecon = newReconBox;
        tempError = tempError + abs(sqrtInt(:,:,j)-abs(detectorWave)).^2;
    end
    errors(iteration) = sum(tempError(:))/totalI;
    dObjectRecon_amp = abs(dObjectRecon);
    dObjectRecon_pha = atan2(imag(dObjectRecon),real(dObjectRecon));
    figure(3),imagesc(dObjectRecon_pha),colorbar;drawnow;
    iteration
    errors(iteration)
    
end




