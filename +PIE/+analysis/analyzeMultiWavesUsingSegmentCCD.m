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
N = 100;% 96, 100
NA = 0.0875;
photon =10000;
stageError_nm = 0;
detSize_mm = 10;
scanningRange_mm =0.0013888;
pie.uieLambda.set(13.56);
pie.uieNA.set(NA);
pie.uieBinning.set(N);
pie.uiez2.set(56.9);
pie.uieDetSize.set(detSize_mm);
pie.uieScanRange.set(scanningRange_mm);
pie.uieGratTilt.set(0);
pie.uieDetTilt.set(0);
pie.uieCenterObstruction.set(0);
pie.uiez1.set(0);
pie.uieNp.set(photon);
pie.uiePhaseShiftingError.set(stageError_nm);
pie.uilSelectMode.setSelectedIndexes(uint8(1));
try
pie.cb(pie.uieLambda);
pie.cb(pie.uieNA);
catch
end
nSteps = 1+round(scanningRange_mm*1e3/pie.do_um(1));
dLambda_nm = [13.12, 13.56, 14.04];
dInt = [0.0099, 1, 0.0113];
dAmp = sqrt(dInt); % probe amp
% dLambda_nm = fliplr(dLambda_nm); % wavelength has to start from short to long 
% dAmp = fliplr(dAmp);

% dAmp=1;
% dLambda_nm =13.56;
modeNumber = length(dLambda_nm);
pie.uieModeNumber.set(modeNumber);

%% path for loading segment 
pie.uieSegmentPath.set(fullfile(pie.cAppPath,  '..','+utils','segs100_2.mat'));

%% load scanning position
scanningDim = '2D';
switch scanningDim
    case '2D'
        % perform 2D scanning
        pie.uieScanSteps.set(nSteps);
        pie.cb(pie.uieScanSteps);
    case '3D'
        % perform 3D scanning
        filename = fullfile(pie.cAppPath,  '..','..', 'data','scanning','squareScanning3D_21by4.mat');
        load(filename);
        pie.dPos_mm = dPos_mm;
        pie.uieScanSteps.set(ceil(sqrt(size(dPos_mm,1))));
        pie.uieScanRange.set(max(max(dPos_mm(:,1:2),[],1)-min(dPos_mm(:,1:2),[],1)));
end


%% probe and object
for u8ModeId = 1:modeNumber
    pie.uilSelectMode.setSelectedIndexes(uint8(u8ModeId));
    % define guessed probe and object
    pie.uieLambda.set(dLambda_nm(u8ModeId));
    pie.uieProbeAmp.set(dAmp(u8ModeId));
    pie.uicbGuess.set(true);
    pie.cb(pie.uibGenProbeObject);
    
    % load zernike aberrations
    % coef = [[4:24]',0.1*randn(21,1)];
    % pie.uieZrn.set(mat2str(coef));
    
    % define probe and object for simulation
    pie.uicbGuess.set(false);
    pie.cb(pie.uibGenProbeObject);
    % load object
    %     pie.cb(pie.uibLoadObject);
    objname = fullfile(pie.cAppPath,  '..','..', 'data','object','contact5_118.mat'); % 5_144, 5_118 contact5_118 threeLine_4pixPitch_145.mat
    load(objname);
    dPosShifts = round((pie.dPos_mm(:,1:2)-min(pie.dPos_mm(:,1:2),[],1))*1000/pie.do_um(u8ModeId));
    K = max(dPosShifts(:,1))+N;
    L = max(dPosShifts(:,2))+N;
    if abs(K-L)<=2
        K=max(K,L);
        L=max(K,L);
    end
    if u8ModeId==1
        pie.dMaxObjectLen = K;
    elseif K<pie.dMaxObjectLen
        K = pie.dMaxObjectLen;
        L = pie.dMaxObjectLen;
    elseif K>pie.dMaxObjectLen
        fprintf('Mode 1 wavelength needs to be the smallest!\n');
        return;
    end
    if K>2000
        fprintf('object sampling: %d, please adjust scanning range\n',K);
        return;
    end
    [m,n] = meshgrid(linspace(0,1,L),linspace(0,1,K));
    [sr,sc]= size(object);
    [p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
    object = interp2(p,q,object,m,n,'nearest');
    try
        pie.dObject(:,:,u8ModeId) = object;
    catch
        pie.dObject = object;
    end
    % Make phase tab active:
    pie.uitgAxesDisplay.selectTabByIndex(pie.U8PROBEOBJECT);
    % Plot wavefronts on phase tab
    pie.replot(pie.U8PROBEOBJECT, []);
    drawnow;
    
end
%% get aerial image
aerialImages_aber = 0;
aerialImages0 = 0;
L = length(pie.dObject(:,:,1));
Lc = 50;
Ls = 1000;
for j =1:length(dLambda_nm)
    [aerialImages,Es] = PIE.utils.getAerialImages(pie.dObject(:,:,j),NA,3000,NA,dLambda_nm(j)/1000,3000,pie.do_um(j),0,L,0);
    aerialImages0 = aerialImages*dInt(j)+aerialImages0;
    pupil_ext = fftshift(fft2(fftshift(pad2(pie.dProbe(:,:,j),L,L))));
    spectrum =  PIE.utils.Propagate (Es,'fourier',pie.dc_um,dLambda_nm(j)/1000,-1);
    spectrum = spectrum.*pupil_ext;
    E_aber =  PIE.utils.Propagate (spectrum,'fourier',pie.dc_um,dLambda_nm(j)/1000,1);
    croppedE = crop2(E_aber,Lc,Lc);
    croppedE = interpFFT(croppedE,Ls,Ls);
    aerialImages_aber = aerialImages_aber+ abs(croppedE).^2;
end

% croppedAerial = interpFFT(croppedAerial,L,L);
xo_um = linspace(-Ls/2,Ls/2,Ls)*pie.do_um(2)/Ls*Lc; % object coordinates
yo_um = linspace(-Ls/2,Ls/2,Ls)*pie.do_um(2)/Ls*Lc; % object coordinates
aerialImages_aber = aerialImages_aber./max(aerialImages_aber(:));
figure(4),imagesc(xo_um,yo_um,aerialImages_aber),xlabel('um'),ylabel('um')
axis equal tight; set(gca,'fontSize',14); colorbar;

%% simulate patterns
pie.cb(pie.uibSimulatePO);

%% reconstruct
pie.uieAlpha.set(0.5);
pie.uieBeta.set(0.03);
pie.uieMaxIteration.set(200);
pie.uieAccuracy.set(0);
pie.uilSelectMode.setSelectedIndexes(uint8(2));
pie.cb(pie.uibComputePhase);

%% analysis
% pie.uilSelectMode.setSelectedIndexes(uint8(2));

pie.uipSelectObject.setSelectedIndex(uint8(3));
pie.uipSelectRegion.setSelectedIndex(uint8(3));
pie.uieSigma.set(0);
pie.cb(pie.uibAnalyze);
resPh = pie.dSelectedObject;
[K,L] = size(resPh);
x_um = pie.dUnit_mm*linspace(-L/2,L/2,L)*1000;
y_um = pie.dUnit_mm*linspace(-K/2,K/2,K)*1000;
resPh(pie.dAnalysisMask==0)=NaN;
resPh =PIE.utils.DelTilt(resPh);
phaseShift = (mean(resPh(pie.dAnalysisMask==1&~isnan(resPh)&angle(pie.dObject(:,:,2))~=0))-...
    mean(resPh(pie.dAnalysisMask==1&~isnan(resPh)&angle(pie.dObject(:,:,2))==0)))/pi
RMS = std(resPh(pie.dAnalysisMask==1&~isnan(resPh)))/pi

axis(pie.haAnalysis,[-inf,inf,-inf,inf,-inf,inf,-1.3,0.7]);
% axis(pie.haAnalysis,[-inf,inf,-inf,inf,-inf,inf,-1,2]);
% figure(2), imagesc(x_um,y_um,resPh);colorbar;axis equal tight;
% xlabel('x/um');ylabel('y/um'); set(gca,'fontSize',14);
L = 1000;
recon = pie.dObjectRecon(:,:,2);
croppedRecon =  crop2(recon,Lc,Lc);
croppedRecon = interpFFT(croppedRecon,L,L);
croppedReconAerial = abs(croppedRecon).^2;
croppedReconAerial = croppedReconAerial./max(croppedReconAerial(:));

figure(5),imagesc(xo_um,yo_um,croppedReconAerial),xlabel('um'),ylabel('um')
axis equal tight; set(gca,'fontSize',14); colorbar;

%% save data
savePath = fullfile(pie.cAppPath,  '..',  '..','data','results',['segment-analysis_', regexprep(datestr(now, 31), ':', '.'),'.mat']);
save(savePath, 'pie','aerialImages0', 'aerialImages_aber','recon');
