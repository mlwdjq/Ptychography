%% this script is used to simulate quadriwave LSI images from EUV mask
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
%% simulation parameters
% basically unchanged parameters
apiPath = 'C:\Program Files\Panoramic\v700';
setupFile = [pie.cAppPath,'\..\..\data\EM-Suite\mask_EUV_QWLSI.sim'];
dataFile = fullfile(pie.cAppPath,'..','..','data','Aerial images',['aerial-images_', regexprep(datestr(now, 31), ':', '.'),'.mat']);
lambda_um = pie.uieLambda.get()/1000; % DUV
Lo_um =pie.uieLo.get()*1000; % object distance
Li_um =pie.uiez2.get()*1000; % imaging distance
scanSteps = pie.uieScanSteps.get(); % scaning steps
scanRange_um = pie.uieScanRange.get(); % scanning range in pupil space
NAo = pie.uieNAo.get(); % Object NA
NAi = sin(atan(Lo_um/Li_um*tan(asin(NAo)))); % imaging NA
det_um = pie.uieDetSize.get()*1000; % detector size
N = length(pie.dObject); % sampling
Nc = pie.uieRes.get(); % sampling
N = Nc;
dc_um = det_um/Nc; % detector size
offsetAngle = 6;
removeTilt =1;
domain_nm = 450; % near field z
t_abs_nm =56; % absorber thickness
method ='TEMPESTpr2';% 'KirchhoffThin';TEMPESTpr2
dz_nm =1;% grid size dz nm
% N = L_nm/dx_nm; % sampling
L_nm = 1000*det_um/(NAo/NAi);
pitch_nm = 400;%L_nm/4/1.882; % feature size
dx_nm =L_nm/N;% grid size dx/dy nm
shearPercentage  = 0.01;
det0_um = 5e3; % real detector size
T_um =  lambda_um./(shearPercentage*det0_um/Li_um/2);
polarDire = 0; % polarization dirction 0 for X-Polarized, 1 for Y-Polarized
 
saveConfig = 1; % save configuration
saveData = 1; % save data
saveImage = 1;
normalIncidence = 0;% normal incident only
filename = 'normalIncidence.mat';
 
% frequently changed parameters
para.domain= domain_nm;
para.lambda = lambda_um*1000; % wavelength nm
para.L = L_nm;
para.pitch = pitch_nm;
para.dx = dx_nm;
para.dz = dz_nm;
para.t_abs = t_abs_nm;
 
%% open simulator
PIE.utils.openSimulator(apiPath);
 
%% load simulation configuration
loadSetup(setupFile);
 
%% change parameters
PIE.utils.setParameters(para,saveConfig,setupFile);
 
%% set illumination
theta_2pi = T_um/2/Lo_um*2;
Nshift = 10;
angleUncertainty = 1e-3;
angleErrorX = angleUncertainty*randn(1,Nshift);
angleErrorY = angleUncertainty*randn(1,Nshift);
deltaTheta = theta_2pi/Nshift;
thetaX = [0:deltaTheta:theta_2pi-deltaTheta]-(theta_2pi-deltaTheta)/2;
thetaY = [0:deltaTheta:theta_2pi-deltaTheta]-(theta_2pi-deltaTheta)/2;
thetaX = thetaX+angleErrorX;
thetaY = thetaY+angleErrorY;
theta = [atand(tand(offsetAngle)+thetaX),atand(sqrt(tand(offsetAngle).^2+thetaY.^2))];
phi = [-atan2d(tand(offsetAngle)+thetaX,0),-atan2d(tand(offsetAngle),-thetaY)];
% phi(phi==0)= 360;
% phi = phi-180;
dA = [theta',phi'];
nInt = length(theta);
if normalIncidence==1
    dA=[6,-90];
    nInt =1;
end
Ex_amp = zeros(N,N,nInt);
Ex_pha = zeros(N,N,nInt);
Ey_amp = zeros(N,N,nInt);
Ey_pha = zeros(N,N,nInt);
aerialImages = cell(nInt,1);
%% run simulator using TEMPEST
for i=1: nInt
    setVariableValues( 'theta',dA(i,1));%
    setVariableValues( 'phi',dA(i,2));%
    
    [~,~,~,Ex_amp(:,:,i),Ex_pha(:,:,i),Ey_amp(:,:,i),Ey_pha(:,:,i)] = PIE.utils.runSimulator(polarDire,method);
    %         pha_near =PIE.utils.UnwrapPhaseBySortingReliabilityWithMask(pha_near,255*ones(N));
end
 
%% get normal incidence
if 1
    setVariableValues( 'theta',offsetAngle);%
    setVariableValues( 'phi',-90);%
    [~,~,~,Ex_ampN,Ex_phaN,Ey_ampN,Ey_phaN] = PIE.utils.runSimulator(polarDire,method);
end
%% get tilt
if removeTilt ==1
    setVariableValues( 'theta',offsetAngle);%
    setVariableValues( 'phi',-90);%
    setVariableValues( 'pitch',0);%
    [~,~,~,Ex_amp0,Ex_pha0,Ey_amp0,Ey_pha0] = PIE.utils.runSimulator(polarDire,method);
    offsetAngle = 0;
else
    Ex_pha0 = 0;
    Ey_pha0 = 0;
end
 
%% generate aerial images
maxInt = 0;
dMaxPhoton = 100;
for i=1: nInt
    switch method
        case 'KirchhoffThin'
            df_um = 0; % defocus distance
        case 'TEMPESTpr2'
            df_um = 0.1; % defocus distance
    end
    Ex = Ex_amp(:,:,i).*exp(1i*Ex_pha(:,:,i))./exp(1i*Ex_pha0);
    Ey = Ey_amp(:,:,i).*exp(1i*Ey_pha(:,:,i))./exp(1i*Ey_pha0);
%     Ey0 = Ey./abs(Ey);
%     Ey0(abs(Ey)<0.6) = 0.5*Ey0(abs(Ey)<0.6);
%     Ey0(abs(Ey)>0.6) = 0.8*Ey0(abs(Ey)>0.6);
    aerialImagesx =  PIE.utils.getQWLSIImages(Ex,NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc,T_um,offsetAngle);
    aerialImagesy =  PIE.utils.getQWLSIImages(Ey,NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc,T_um,offsetAngle);
    aerialImages{i}= aerialImagesx + aerialImagesy;
    
    % photon noise
    wInt = aerialImages{i};
    if (dMaxPhoton > 0)
        wInt = wInt * dMaxPhoton;
        wInt = wInt + sqrt(wInt).*randn(size(wInt));
        wInt = wInt / dMaxPhoton;
        aerialImages{i} = wInt;
    end
    
       imagesc(aerialImages{i}),axis xy; colorbar;pause(0.5);
%             abs(Ey(100,100))
    if max(max(abs(aerialImages{i})))>maxInt
        maxInt = max(max(aerialImages{i}));
    end
end
 
%% get origin phase
EN = Ey_ampN.*exp(1i*Ey_phaN)./exp(1i*Ey_pha0);
[~,EN] =  PIE.utils.getAerialImages(EN,NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc,offsetAngle);
Ex_phaN = -unwrap(unwrap(angle(EN),[],1),[],2)/2/pi;
% %% save normal incidence object
% if normalIncidence==1
%     object = PIE.utils.Propagate (Ey,'angular spectrum',dx_nm/1000,lambda_um,df_um);
%     save([pie.cAppPath,'\..\..\data\object\',filename],'object');
%     return;
% end
 
%% save data
if saveData==1
    save(dataFile, 'para','Ex_amp','Ey_amp','Ex_pha','Ey_pha',...
        'Ex_amp0','Ex_pha0','Ey_amp0','Ey_pha0', 'Ex_ampN','Ex_phaN',...
        'Ey_ampN','Ey_phaN','theta', 'phi','aerialImages','angleErrorX','angleErrorY');
end
%% save images
if saveImage == 1
    mkdir(dataFile(1:end-4));
    for i= 1:nInt
         imwrite(aerialImages{i}/maxInt,fullfile(dataFile(1:end-4),[num2str(i),'.bmp']));
%         aerialImagesScale{i} = aerialImages{i}/maxInt;
%         aerialImagesScale{i} =aerialImages{i};
%          imagesc(aerialImagesScale{i}-aerialImagesMo{i}),axis xy; colorbar;pause(0.5);
%          imwrite(mat2gray(aerialImagesScale{i}-aerialImagesMo{i}),[dataFile(1:end-4),'\',num2str(i),'.bmp']);
    end
end
 
%% reconstruct
PIE.simulation.reconstructQWLSI();
 
%% plot result
% figure(2),imagesc(aerialImages(:,:,1)), axis off tight equal; colorbar;
