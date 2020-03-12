%% this script is used to simulate aerial images from EUV mask

%% simulation parameters
% basically unchanged parameters
apiPath = 'C:\Program Files\Panoramic\v700';
setupFile = 'D:\OneDrive\Ptychography\code\Ptychography\data\EM-Suite\mask_EUV2.sim';
dataFile = ['D:\OneDrive\Ptychography\code\Ptychography\data\Aerial images\aerial-images_', regexprep(datestr(now, 31), ':', '.'),'.mat'];
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
df_um =0; % defocus distance
domain_nm = 450; % near field z 
t_abs_nm =60; % absorber thickness
method ='TEMPESTpr2';% 'KirchhoffThin';TEMPESTpr2
dz_nm =1;% grid size dz nm
% N = L_nm/dx_nm; % sampling
L_nm = 1000*det_um/(NAo/NAi);
pitch_nm = L_nm/4; % feature size
dx_nm =L_nm/N;% grid size dx/dy nm

polarDire = 0; % polarization dirction 0 for X-Polarized, 1 for Y-Polarized

saveConfig = 1; % save configuration
saveData = 1; % save data
saveImage = 1;


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
load('D:\OneDrive\Ptychography\code\Ptychography\data\scanning\circular_FPM_seg16.mat');
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
E = zeros(N,N,scanSteps^2);
aerialImages = cell(scanSteps^2,1);
maxInt = 0;
for i=1: nInt
        setVariableValues( 'theta',dA(i,1));
        setVariableValues( 'phi',dA(i,2));
        %% run simulator using TEMPEST
        E(:,:,i) = PIE.utils.runSimulator(polarDire,method);
        %         pha_near =PIE.utils.UnwrapPhaseBySortingReliabilityWithMask(pha_near,255*ones(N));
        if strcmp(method,'KirchhoffThin')
            E(:,:,i) = E(:,:,i)./abs(E(1,1,i));
        end
        %% generate aerial images
        aerialImages{i}=  PIE.utils.getAerialImages(E(:,:,i),NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc);
        if max(max(abs(aerialImages{i})))>maxInt
            maxInt = max(max(aerialImages{i}));
        end
end
%% remove defocus
% df_um =0;
% maxInt=0;
% for i=1: nInt
%     aerialImages{i}=  PIE.utils.getAerialImages(E(:,:,i),NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc);
%     figure(2),imagesc(aerialImages{i});pause(0.1);
%     if max(max(abs(aerialImages{i})))>maxInt
%         maxInt = max(max(aerialImages{i}));
%     end
% end

%% save data
if saveData==1
    save(dataFile, 'para', 'E', 'theta', 'phi','aerialImages','dPos_mm');
end
%% save images
if saveImage == 1
    mkdir(dataFile(1:end-4));
    for i= 1:nInt
        imwrite(aerialImages{i}/maxInt,[dataFile(1:end-4),'\',num2str(i),'.bmp']);
    end
end

%% plot result
% figure(2),imagesc(aerialImages(:,:,1)), axis off tight equal; colorbar;