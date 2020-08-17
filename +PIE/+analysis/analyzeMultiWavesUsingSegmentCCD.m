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
NA = 0.0875;
pie.uieLambda.set(13.5);
pie.uieNA.set(NA);
pie.uieBinning.set(N);
pie.uiez2.set(17);
pie.uieDetSize.set(10);
pie.uieScanRange.set(0.000459);
pie.uieBinning.set(96);
pie.uieGratTilt.set(0);
pie.uieDetTilt.set(0);
pie.uieCenterObstruction.set(0);
pie.uiez1.set(0);
pie.uieNp.set(0);

dLambda_nm = [16.33, 15.69, 15.09, 14.55, 14.04, 13.56,...
    13.11, 12.07, 12.31, 11.94, 11.59];

dInt = [6.08E-02, 1.34E-01, 2.44E-01, 3.95E-01, 6.50E-01, 1.00E+00,...
    1.30E+00, 1.11E+00, 4.22E-01, 1.55E-01, 3.32E-02];
dAmp = sqrt(dInt); % probe amp
dLambda_nm = fliplr(dLambda_nm); % wavelength has to start from short to long 
dAmp = fliplr(dAmp);

% dAmp=1;
% dLambda_nm =13.5;
modeNumber = length(dLambda_nm);
pie.uieModeNumber.set(modeNumber);

%% path for loading segment 
pie.uieSegmentPath.set(fullfile(pie.cAppPath,  '..','+utils','segs32.mat'));

%% load scanning position
scanningDim = '2D';
switch scanningDim
    case '2D'
        % perform 2D scanning
        pie.uieScanSteps.set(21);
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
    objname = fullfile(pie.cAppPath,  '..','..', 'data','object','contact11_115.mat');
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
    pie.dObject(:,:,u8ModeId) = object;
    % Make phase tab active:
    pie.uitgAxesDisplay.selectTabByIndex(pie.U8PROBEOBJECT);
    % Plot wavefronts on phase tab
    pie.replot(pie.U8PROBEOBJECT, []);
    drawnow;
    
end
%% simulate patterns
pie.cb(pie.uibSimulatePO);

%% reconstruct
pie.uieAlpha.set(0.2);
pie.uieBeta.set(0.02);
pie.uieMaxIteration.set(200);
pie.uieAccuracy.set(0);
pie.cb(pie.uibComputePhase);

%% analysis
pie.uilSelectMode.setSelectedIndexes(uint8(6));

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
RMS = std(resPh(pie.dAnalysisMask==1&~isnan(resPh)))
figure(2), imagesc(x_um,y_um,resPh);colorbar;axis equal tight;
xlabel('x/um');ylabel('y/um'); set(gca,'fontSize',14);

