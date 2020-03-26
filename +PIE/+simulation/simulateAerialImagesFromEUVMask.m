%% this script is used to simulate aerial images from EUV mask

%% simulation parameters
% basically unchanged parameters
apiPath = 'C:\Program Files\Panoramic\v700';
setupFile = [pie.cAppPath,'\..\..\data\EM-Suite\mask_EUV2.sim'];
dataFile = [pie.cAppPath,'\..\..\data\Aerial images\aerial-images_', regexprep(datestr(now, 31), ':', '.'),'.mat'];
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
dc_um = det_um/Nc; % detector size
offsetAngle = 6.1605;
domain_nm = 450; % near field z
t_abs_nm =56; % absorber thickness
method ='KirchhoffThin';% 'KirchhoffThin';TEMPESTpr2
dz_nm =1;% grid size dz nm
% N = L_nm/dx_nm; % sampling
L_nm = 1000*det_um/(NAo/NAi);
pitch_nm = L_nm/4; % feature size
dx_nm =L_nm/N;% grid size dx/dy nm

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
dPos = -scanRange_um/2:scanRange_um/(scanSteps-1):scanRange_um/2;
[dm,dn] = meshgrid(dPos,dPos);
dPos_um  = [dm(:),dn(:)];
dPos_mm = dPos_um/1000;
load([pie.cAppPath,'\..\..\data\scanning\CircularRoundScanningSeg16.mat']);
nInt = size(dPos_mm,1);
theta= atan(sqrt(dPos_mm(:,1).^2+dPos_mm(:,2).^2)/(Lo_um/1000))/pi*180; % illumination angles
phi = atan2(-dPos_mm(:,2),-dPos_mm(:,1))/pi*180; % illumination azimuthes
dA = [theta(:),phi(:)];
% ky = -sin(atan(dPos_mm(:,1)/(Lo_um/1000)));
% kx = -sin(atan(dPos_mm(:,2)/(Lo_um/1000)));
% kx2 = -sin(theta/180*pi).*sin(phi/180*pi);
% ky2 = -sin(theta/180*pi).*cos(phi/180*pi);
% [x,y] = meshgrid(linspace(-1,1,150));
% pha = 2*pi*kx2(3).*x*25.49+2*pi*ky2(3).*y*25.49;
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


%% generate aerial images
maxInt = 0;
for i=1: nInt
    switch method
        case 'KirchhoffThin'
            df_um = 0; % defocus distance
            %                 Ex = Ex_amp(:,:,i).*exp(1i*Ex_pha(:,:,i));
            %                 total = sum(abs(Ex(:)).^2);
            %                 Ex = Ex/sqrt(total);
            %                 aerialImages{i}=  PIE.utils.getAerialImages(Ex,NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc,offsetAngle);
            Ex = Ex_amp(:,:,i).*exp(1i*Ex_pha(:,:,i));
            Ey = Ey_amp(:,:,i).*exp(1i*Ey_pha(:,:,i));
            %                 total = sum(abs(Ex(:)).^2)+sum(abs(Ey(:)).^2);
            %                 Ex = Ex/sqrt(total);
            %                 Ey = Ey/sqrt(total);
            aerialImagesx =  PIE.utils.getAerialImages(Ex,NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc,offsetAngle);
            aerialImagesy =  PIE.utils.getAerialImages(Ey,NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc,offsetAngle);
            aerialImages{i}= aerialImagesx+aerialImagesy;
        case 'TEMPESTpr2'
            df_um = 0.1; % defocus distance
            Ex = Ex_amp(:,:,i).*exp(1i*Ex_pha(:,:,i));
            Ey = Ey_amp(:,:,i).*exp(1i*Ey_pha(:,:,i));
            %                 total = sum(abs(Ex(:)).^2)+sum(abs(Ey(:)).^2);
            %                 Ex = Ex/sqrt(total);
            %                 Ey = Ey/sqrt(total);
            aerialImagesx =  PIE.utils.getAerialImages(Ex,NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc,offsetAngle);
            aerialImagesy =  PIE.utils.getAerialImages(Ey,NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc,offsetAngle);
            aerialImages{i}= aerialImagesx+aerialImagesy;
    end
%     imagesc(aerialImages{i}),axis xy; colorbar;pause(0.5);
    if max(max(abs(aerialImages{i})))>maxInt
        maxInt = max(max(aerialImages{i}));
    end
end

%% save normal incidence object
if normalIncidence==1
    object = PIE.utils.Propagate (Ey,'angular spectrum',dx_nm/1000,lambda_um,df_um);
    save([pie.cAppPath,'\..\..\data\object\',filename],'object');
    return;
end

%% save data
if saveData==1
    save(dataFile, 'para','Ex_amp','Ey_amp','Ex_pha','Ey_pha', 'theta', 'phi','aerialImages','dPos_mm');
end
%% save images
if saveImage == 1
    mkdir(dataFile(1:end-4));
    for i= 1:nInt
         imwrite(aerialImages{i}/maxInt,[dataFile(1:end-4),'\',num2str(i),'.bmp']);
%         aerialImagesScale{i} = aerialImages{i}/maxInt;
%         aerialImagesScale{i} =aerialImages{i};
%          imagesc(aerialImagesScale{i}-aerialImagesMo{i}),axis xy; colorbar;pause(0.5);
%          imwrite(mat2gray(aerialImagesScale{i}-aerialImagesMo{i}),[dataFile(1:end-4),'\',num2str(i),'.bmp']);
    end
end

%% plot result
% figure(2),imagesc(aerialImages(:,:,1)), axis off tight equal; colorbar;