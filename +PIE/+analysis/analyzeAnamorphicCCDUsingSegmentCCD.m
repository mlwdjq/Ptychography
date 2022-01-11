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
N = 96;
NA = 0.55/4;
NAy = 0.55/8;
photon =10000;
stageError_nm = 0;
detSize_mm = 40;
scanningRange_mm =0.0017418;
pie.uieLambda.set(13.56);
pie.uieNA.set(NA);
pie.uieProbeEllipticity.set(NAy/NA);
pie.uieBinning.set(N);
pie.uiez2.set(73.4);
pie.uieDetSize.set(detSize_mm);
pie.uieScanRange.set(scanningRange_mm);
pie.uieGratTilt.set(0);
pie.uieDetTilt.set(0);
pie.uieCenterObstruction.set(0);
pie.uiez1.set(0);
pie.uieNp.set(photon);
pie.uiePhaseShiftingError.set(stageError_nm);
pie.uilSelectMode.setSelectedIndexes(uint8(1));
pie.uipSelectMask.setSelectedIndex(uint8(1));
pie.cb(pie.uipSelectMask);% anamorphic CCD 2:1
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
pie.uieSegmentPath.set(fullfile(pie.cAppPath,  '..','+utils','segs16e_2.mat'));

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
    objname = fullfile(pie.cAppPath,  '..','..', 'data','object','anamorphicLine_136_2.mat');
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


%% simulate patterns
pie.cb(pie.uibSimulatePO);

%% reconstruct
pie.uieAlpha.set(0.5);
pie.uieBeta.set(0.02);
pie.uieMaxIteration.set(200);
pie.uieAccuracy.set(0);
pie.uilSelectMode.setSelectedIndexes(uint8(2));
pie.cb(pie.uibComputePhase);

%% analysis
pie.uipSelectObject.setSelectedIndex(uint8(12));
pie.uipSelectRegion.setSelectedIndex(uint8(3));
pie.uieSigma.set(3);
pie.cb(pie.uibAnalyze);
resPh = pie.dSelectedObject;
[K,L] = size(resPh);
x_um = pie.dUnit_mm*linspace(-L/2,L/2,L)*1000;
y_um = pie.dUnit_mm*linspace(-K/2,K/2,K)*1000;
resPh(pie.dAnalysisMask==0)=NaN;
resPh =PIE.utils.DelTilt(resPh);
RMS = std(resPh(pie.dAnalysisMask==1&~isnan(resPh)))/pi
% figure(2), imagesc(x_um,y_um,resPh);colorbar;axis equal tight;
% xlabel('x/um');ylabel('y/um'); set(gca,'fontSize',14);
%%
dObjectRecon = pie.dObjectRecon(:,:,2);
dObject = pie.dObject(:,:,2);
N = length(dObject);
dObject_upSampling{1} = pad2(dObject,2*N,N);
dObject_upSampling{2} = pad2(dObjectRecon,2*N,N);

for i = 1:2
    % generate defocused aerial image
    dObjectScale = dObject_upSampling{i};
    df_um = 0.5;
    z_um       = pie.uiez2.get()*1000;
    Rc_um = (z_um+df_um)*tan(asin(NA));
    [n1,n2]=meshgrid(1:N);
    n1 = n1-N/2-1;
    n2 = n2-N/2-1;
    lambda_um = dLambda_nm(2)/1000;
    dc_um = lambda_um*z_um/N/pie.do_um(2);
    Hs= exp(-1i*pi*df_um*dc_um^2/lambda_um/z_um^2*(n1.^2+n2.^2));% defocus
    pupil = pinhole(round(2*Rc_um/dc_um),N,N).*Hs;
    spectrum =  PIE.utils.Propagate (dObjectScale,'fourier',pie.do_um(2),lambda_um,-1);
    spectrum = crop2(spectrum,N,N).*pupil;%imagesc(abs(spectrum));
    Ns = 10*N; % increase sampling
    spectrum = pad2(spectrum,Ns,Ns);
    Es =  PIE.utils.Propagate (spectrum,'fourier',pie.do_um(2),lambda_um,1);
    aerialImages = abs(Es).^2;
    aerialImages = mat2gray(aerialImages);
    aerialImagess{i} = aerialImages;
    % crop for a better view
    Nc = Ns/4;
    x_um = linspace(-N/2,N/2,Nc)*pie.do_um(2)/4;
    figure(i+5),imagesc(x_um,x_um,crop2(aerialImages,Nc,Nc));axis equal tight; xlabel('x/um');ylabel('y/um');colorbar;
    
end
% direction imaging
dObjectScale = dObjectRecon;
spectrum =  PIE.utils.Propagate (dObjectScale,'fourier',pie.do_um(2),lambda_um,-1);
pupil2 = lsianalyze.utils.elipticalHole(round(2*Rc_um/dc_um),...
    round(2*Rc_um/dc_um*tan(asin(NAy))/tan(asin(NA))),N,N).*Hs;
spectrum = spectrum.*pupil2;%imagesc(abs(spectrum));
Ns = 10*N; % increase sampling
spectrum = pad2(spectrum,Ns/2,Ns);
Es =  PIE.utils.Propagate (spectrum,'fourier',pie.do_um(2),lambda_um,1);
Es = pad2(Es,Ns,Ns);
aerialImages = abs(Es).^2;
aerialImages = mat2gray(aerialImages);
aerialImagess{3} = aerialImages;
% crop for a better view
Nc = Ns/4;
x_um = linspace(-N/2,N/2,Nc)*pie.do_um(2)/4;
figure(4),imagesc(x_um,x_um,crop2(aerialImages,Nc,Nc));axis equal tight; xlabel('x/um');ylabel('y/um');colorbar;

% figure(5),imagesc(x_um,x_um,crop2(aerialImagess{2}-aerialImagess{1},Nc,Nc));axis equal tight; xlabel('x/um');ylabel('y/um');colorbar;