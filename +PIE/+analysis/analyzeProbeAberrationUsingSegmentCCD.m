%% this script is used to evaluate the probe reconstruction performance from PIE
try
    if ~exist('pie')||~ishandle(pie.hFigure)
        launch_PIE;
        msgbox('Please change PIE paras setup and run again');
        return;
    end
catch
    launch_PIE;
    msgbox('Please change PIE paras setup and run again');
    return;
end

%% setup
% remember loading a setup on the GUI

%% path for loading segment 
pie.uieSegmentPath.set(fullfile(pie.cAppPath,  '..','+utils','segs32.mat'));

%% load scanning position
filename = fullfile(pie.cAppPath,  '..','..', 'data','scanning','CircularScanningSeg_N9Nz11.mat');
load(filename);
pie.dPos_mm = dPos_mm;
pie.uieScanSteps.set(ceil(sqrt(size(dPos_mm,1))));
pie.uieScanRange.set(max(max(dPos_mm(:,1:2),[],1)-min(dPos_mm(:,1:2),[],1)));

% perform 2D scanning
% pie.uieScanSteps.set(16);
% pie.cb(pie.uieScanSteps);

%% define guessed probe and object
pie.uieZrn.set('[]');
pie.uicbGuess.set(true);
pie.cb(pie.uibGenProbeObject);

%% load zernike aberrations
coef = [[4:24]',0.1*randn(21,1)];
pie.uieZrn.set(mat2str(coef));

%% define probe and object for simulation
pie.uicbGuess.set(false);
pie.cb(pie.uibGenProbeObject);
pie.cb(pie.uibLoadObject);

%% simulate patterns
pie.cb(pie.uibSimulatePO);

%% reconstruct
pie.uieAlpha.set(0.5);
pie.uieBeta.set(0.2);
pie.uieMaxIteration.set(100);
pie.uieAccuracy.set(0);
pie.cb(pie.uibComputePhase);

%% analysis
Rc_um   = pie.uieRprobe.get()*1000;
N           = pie.uieRes.get();
mask = pinhole(round(2*Rc_um/pie.dc_um),N,N);
pupil = fftshift(fft2(fftshift(pie.dProbeRecon)));
pupil0 = fftshift(fft2(fftshift(pie.dProbe)));
phase = angle(pupil);
phase = phase.*mask/2/pi;
phase0 = angle(pupil0);
phase0 = phase0.*mask/2/pi;
res = (phase-phase0);
res = res - mean(res(mask==1));
rms = std(res(mask==1));
figure(4),imagesc(phase);colorbar; axis tight equal off;set(gca,'fontSize',14);
figure(5),imagesc(phase0);colorbar; axis tight equal off;set(gca,'fontSize',14);
figure(6),imagesc(res);colorbar; axis tight equal off;set(gca,'fontSize',14);

% generate zernike basis
Nz =24;
detSize_um  = pie.uieDetSize.get()*1000;
[x_um,y_um] = meshgrid(linspace(-detSize_um/2,detSize_um/2,N));
df_um       = pie.uiez1.get()*1000;% negative sign corresponds convergent
z_um       = pie.uiez2.get()*1000;
NA           = pie.uieNA.get();
r = sin(atan(sqrt(x_um.^2+y_um.^2)/(z_um+df_um)))/NA;
th = atan2(y_um,x_um);
rs = r(mask ==1);
ths = th(mask==1);
basis = zeros(length(rs),Nz+1);
for k =0:Nz
    afn = zgen([], k, 'fnr');
    basis(:,k+1)=afn(rs, ths);
end
dZrn0 = pinv(basis)*phase0(mask==1)/2/pi;
dZrn = pinv(basis)*phase(mask==1)/2/pi;
dZrnRes =dZrn-dZrn0;
dZrnRes(1:4) =0;

figure(3),h =plot(0:Nz,dZrn);xlabel('Zernike terms');ylabel('Residual errors/waves');
set(h,'lineWidth',2);set(gca,'fontSize',14);hold on;
figure(3),h =plot(0:Nz,dZrn0);xlabel('Zernike terms');ylabel('Residual errors/waves');
set(h,'lineWidth',2);set(gca,'fontSize',14);hold on;
figure(3),h =plot(0:Nz,dZrnRes);xlabel('Zernike terms');ylabel('Residual errors/waves');
set(h,'lineWidth',2);set(gca,'fontSize',14);hold on;

