%% this script is used to simulate aerial images from EUV mask
clc;
clear;

%% simulation parameters
% basically unchanged parameters
apiPath = 'C:\Program Files\Panoramic\v700';
setupFile = 'D:\OneDrive\Ptychography\code\Ptychography\data\EM-Suite\mask_EUV.sim';
dataFile = ['D:\OneDrive\Ptychography\code\Ptychography\data\Aerial images\aerial-images_', regexprep(datestr(now, 31), ':', '.'),'.mat'];
lambda_um = 13.5e-3; % DUV
Lo_um =3000; % object distance
Li_um =1000000; % imaging distance
scanSteps = 4; % scaning steps
scanRange_um = 2000; % scanning range in pupil space
NAo = 0.5; % Object NA
NAi = sin(atan(Lo_um/Li_um*tan(asin(NAo)))); % imaging NA
nD = 50; % far field sampling
det_um = 200; % detector size
N = 150; % sampling
dc_um = det_um/N; % detector size
df_um =0; % defocus distance
pitch_nm = 150; % feature size
domian_nm = 450; % near field z 
dx_nm =4;% grid size dx/dy nm
dz_nm =1;% grid size dz nm
t_Si_nm = 3.05; % Si thickness
t_Mo_nm = 2.4; % Mo thickness
% N = L_nm/dx_nm; % sampling
L_nm = 1000*det_um/(Li_um/Lo_um);

polarDire = 0; % polarization dirction 0 for X-Polarized, 1 for Y-Polarized

saveConfig = 1; % save configuration
saveData = 1; % save data


% frequently changed parameters
para.domian = domian_nm; 
para.lambda = lambda_um*1000; % wavelength nm
para.L = L_nm;
para.pitch = pitch_nm; 
para.dx = dx_nm; 
para.dz = dz_nm; 
para.t_Si=t_Si_nm;
para.t_Mo=t_Mo_nm;

%% open simulator
PIE.utils.openSimulator(apiPath);

%% load simulation configuration
loadSetup(setupFile);

%% change parameters
PIE.utils.setParameters(para,saveConfig,setupFile);

%% set illumination
dPos = 0:scanRange_um/(scanSteps-1):scanRange_um;
[dm,dn] = meshgrid(dPos,dPos);
dPos_um  = [dm(:),dn(:)];
dPos_mm = dPos_um/1000;
[dm,dn] = meshgrid(dPos-scanRange_um/2,dPos-scanRange_um/2);
theta= atan(sqrt(dm.^2+dn.^2)/Lo_um)/pi*180; % illumination angles
phi = atan2(dm,dn)/pi*180; % illumination azimuthes
dA = [theta(:),phi(:)];

E = zeros(N,N,scanSteps^2);
aerialImages = cell(scanSteps^2,1);
for i=1: scanSteps^2
        setVariableValues( 'theta',dA(i,1));
        setVariableValues( 'phi',dA(i,2));
        %% run simulator using TEMPEST
        E(:,:,i) = PIE.utils.runSimulator(polarDire);
        %% generate aerial images
        aerialImages{i}=  PIE.utils.getAerialImages(E(:,:,i),Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,N);
end

%% save data
if saveData==1
    save(dataFile, 'para', 'E', 'theta', 'phi','aerialImages','dPos_mm');
    for i= 1:scanSteps^2
       imwrite(mat2gray(aerialImages{i}),[dataFile(1:end-1),'_',num2str(i),'.bmp']);
    end
end

%% plot result
% figure(2),imagesc(aerialImages(:,:,1)), axis off tight equal; colorbar;