%% this script is used to analyze the incident aberration effects of QWLSI

%% load data
% load('aerial-images_2020-07-24 09.02.45.mat');
lambda_um = pie.uieLambda.get()/1000; % DUV
Lo_um =pie.uieLo.get()*1000; % object distance
Li_um =pie.uiez2.get()*1000; % imaging distance
NAo = pie.uieNAo.get(); % Object NA
NAi = sin(atan(Lo_um/Li_um*tan(asin(NAo)))); % imaging NA
det_um = pie.uieDetSize.get()*1000; % detector size
Nc = pie.uieRes.get(); % sampling
dc_um = det_um/Nc; % detector size
shearPercentage  = 0.01;
det0_um = 5e3; % real detector size
L_nm = 1000*det_um/(NAo/NAi);
T_um =  lambda_um./(shearPercentage*det0_um/Li_um/2);
%% define aberrations
try 
    noAber;
catch
    [N,~,nInt] = size(Ex_amp);
    Nz = 8;
    [x,y] = meshgrid(linspace(-1.5,1.5,N));
    [th,r] = cart2pol(x,y);
    for k = 3:Nz
        afn = zgen([], k, 'fnr');
        temp=afn(r, th);
        Basis(:,k-2) = temp(:);
    end
    aber = zeros(N,N,Nz);
    
    for i = 1:nInt
        dZrn = randn(Nz-2,1)*0.1;
        temp = reshape(Basis*dZrn,N,N);
        aber(:,:,i) = temp;
    end
end
%% generate aerial images
method = 'TEMPESTpr2';
dMaxPhoton =0;
for i=1: nInt
    switch method
        case 'KirchhoffThin'
            df_um = 0; % defocus distance
        case 'TEMPESTpr2'
            df_um = 0.1; % defocus distance
    end
    Ex = Ex_amp(:,:,i).*exp(1i*Ex_pha(:,:,i))./exp(1i*Ex_pha0).*exp(1i*2*pi*aber(:,:,i));
    Ey = Ey_amp(:,:,i).*exp(1i*Ey_pha(:,:,i))./exp(1i*Ey_pha0).*exp(1i*2*pi*aber(:,:,i));
    %     Ey0 = Ey./abs(Ey);
    %     Ey0(abs(Ey)<0.6) = 0.5*Ey0(abs(Ey)<0.6);
    %     Ey0(abs(Ey)>0.6) = 0.8*Ey0(abs(Ey)>0.6);
    aerialImagesx =  PIE.utils.getQWLSIImages(Ex,NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc,T_um,0);
    aerialImagesy =  PIE.utils.getQWLSIImages(Ey,NAo,Lo_um,NAi,lambda_um,Li_um,dc_um,df_um,Nc,T_um,0);
    aerialImages{i}= aerialImagesx + aerialImagesy;
    % photon noise
    wInt = aerialImages{i};
    if (dMaxPhoton > 0)
        wInt = wInt * dMaxPhoton;
        wInt = wInt + sqrt(wInt).*randn(size(wInt));
        wInt = wInt / dMaxPhoton;
        aerialImages{i} = wInt;
    end
    %              imagesc(aerialImages{i}),axis xy; colorbar;pause(0.5);
end

%% reconstruct
PIE.simulation.reconstructQWLSI();