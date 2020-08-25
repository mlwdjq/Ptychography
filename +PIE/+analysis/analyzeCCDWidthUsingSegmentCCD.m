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
N = 100;
NA = 0.0875;
photon =10000;
detSize_mm = 30;
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
pie.uilSelectMode.setSelectedIndexes(uint8(1));
pie.cb(pie.uieLambda);
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
    objname = fullfile(pie.cAppPath,  '..','..', 'data','object','contact15_155.mat');
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
pie.uieBeta.set(0.03);
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

