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
N = 128;
NA = 0.00026351;
lambda_um =13.5e-3;
pie.uicbFourierPtychography.set(true);
pie.uieLambda.set(lambda_um*1e3);
pie.uieNA.set(NA);
pie.uieBinning.set(N);
pie.uiez2.set(1000);
pie.uieDetSize.set(0.2);% 0.2, 0.5s
pie.uieScanRange.set(0.4282);
pie.uieBinning.set(N);
pie.uieGratTilt.set(0);
pie.uieDetTilt.set(0);
pie.uieCenterObstruction.set(0);
pie.uiez1.set(0);
pie.uieNp.set(100);% 100, 400
pie.cb(pie.uieNA);

modeNumber = 1;
pie.uieModeNumber.set(modeNumber);


%% load scanning position

filename = fullfile(pie.cAppPath,  '..','..', 'data','scanning','CircularRoundOffsetScanningSeg16s-6.4188.mat');%6.16405, 6.4188
load(filename);
pie.dPos_mm = dPos_mm;
pie.uieScanSteps.set(ceil(sqrt(size(dPos_mm,1))));
pie.uieScanRange.set(max(max(dPos_mm(:,1:2),[],1)-min(dPos_mm(:,1:2),[],1)));



%% probe and object
for u8ModeId = 1:modeNumber
    pie.uilSelectMode.setSelectedIndexes(uint8(u8ModeId));
    % define guessed probe and object
%     pie.uieLambda.set(dLambda_nm(u8ModeId));
%     pie.uieProbeAmp.set(dAmp(u8ModeId));
    pie.uieProbeOffset.set('[-0.3375,0]');%'[-0.324,0]' '[-0.3375,0]'
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
    
    objname = fullfile(pie.cAppPath,  '..','..', 'data','object','56nm-Mo-6.4188.mat');
    load(objname);
    
    k0 = 2*pi/lambda_um;
    kr = k0*sin(atan(sqrt(pie.dPos_mm(:,1).^2+pie.dPos_mm(:,2).^2)/pie.uieLo.get()));
    phi = atan2(pie.dPos_mm(:,1),pie.dPos_mm(:,2));
    kxy = [kr.*sin(phi),kr.*cos(phi)];
    dkxy = 2*pi/pie.dc_um/N*pie.uieMag.get();
    dPosShifts = kxy./dkxy;
    dPosShifts = dPosShifts -min(dPosShifts,[],1);
    dPosShifts = round(dPosShifts);
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
    object =  PIE.utils.Propagate (object,'fourier',pie.do_um(u8ModeId),lambda_um,-1);
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
pie.uieAlpha.set(0.5);
pie.uieBeta.set(0.05);
pie.uieMaxIteration.set(200);
pie.uieAccuracy.set(0);
pie.cb(pie.uibComputePhase);

%% analysis
CTF = zeros(K,L);
for i =1:length(dPosShifts)
    CTF(dPosShifts(i,1)+[1:N],dPosShifts(i,2)+[1:N]) = pie.dProbe + CTF(dPosShifts(i,1)+[1:N],dPosShifts(i,2)+[1:N]);
end
CTF(CTF~=0)=1;
% pie.dObject = pad2(crop2(object,N,N).*pie.dProbe,K,L);
pie.dObject =object.*CTF;
pie.uilSelectMode.setSelectedIndexes(uint8(1));

pie.uipSelectObject.setSelectedIndex(uint8(12));
pie.uipSelectRegion.setSelectedIndex(uint8(1));
pie.uieSigma.set(3);
pie.cb(pie.uibAnalyze);
resPh = pie.dSelectedObject;
[K,L] = size(resPh);
x_um = pie.dUnit_mm*linspace(-L/2,L/2,L)*1000;
y_um = pie.dUnit_mm*linspace(-K/2,K/2,K)*1000;
resPh(pie.dAnalysisMask==0)=NaN;
resPh =PIE.utils.DelTilt(resPh);
RMS = std(resPh(pie.dAnalysisMask==1&~isnan(resPh)))/pi
figure(2), imagesc(x_um,y_um,resPh);colorbar;axis equal tight xy;
xlabel('x/um');ylabel('y/um'); set(gca,'fontSize',14);hold on;
mask2 = zeros(size(resPh));
rh =5;
mask2(N/2+1+rh:N/4*3-rh,N/2+1+rh:N/4*3-rh)=1;
avg = mean(resPh(mask2==1&~isnan(resPh)))/pi
% imagesc(x_um,y_um,mask2);

%% save aerial images
saveImage =0;
dataFile = [pie.cAppPath,'\..\..\data\Aerial images\aerial-images_', regexprep(datestr(now, 31), ':', '.'),'.mat'];
if saveImage == 1
    mkdir(dataFile(1:end-4));
    maxInt=0;
    for i = 1:length(dPosShifts)
        if (max(pie.ceInt{i}(:))>maxInt)
            maxInt = max(pie.ceInt{i}(:));
        end
    end
    for i= 1:length(dPosShifts)
         imwrite(pie.ceInt{i}/maxInt,[dataFile(1:end-4),'\',num2str(i),'.bmp']);
    end
end