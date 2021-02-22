%% this script is used to analyze the incident aberration effects of QWLSI by calibrate the null test result


%% load data
load('aerial-images_2021-02-21 23.12.39 (QWLSI_s0.75um_null).mat');
[N,~,nInt] = size(Ex_amp);

%% define aberration
noAber = 0;
aberAmp = linspace(0,0.1,11);
% aberAmp = 0.1;
RMS = zeros(1, length(aberAmp));
dWs_PS = zeros(1, length(aberAmp));
Nz = 16;
[x,y] = meshgrid(linspace(-0.7,0.7,N));
[th,r] = cart2pol(x,y);
for k =5:Nz
    afn = zgen([], k, 'fnr');
    temp=afn(r, th);
    Basis(:,k) = temp(:);
end
aber = zeros(N,N,Nz);

for i = 1:nInt
    dZrn = randn(Nz,1);
    temp = reshape(Basis*dZrn,N,N);
    aber(:,:,i) = temp./std(temp(:));
%     figure(3),imagesc(aber(:,:,i));pause(1)
end
aber0 = aber;
for j = 1:length(aberAmp)
    %% load data
    load('aerial-images_2021-02-21 23.12.39 (QWLSI_s0.75um_null).mat');
%     [N,~,nInt] = size(Ex_amp);
aberAmp(j)
    aber = aber0*aberAmp(j);
    
    % analysis
    PIE.analysis.QWLSIIncidentAberrationAnalysis;
    dWx0 = dWx;
    dWy0 = dWy;
    
    %% load null test data
    load('aerial-images_2021-02-21 23.06.43 (QWLSI_s0.75um).mat');
    PIE.analysis.QWLSIIncidentAberrationAnalysis;
    dWx1 = dWx;
    dWy1 = dWy;
    
    %% reconstruction
    dWxss = exp(1i*dWx1)./exp(1i*dWx0);
    dWyss = exp(1i*dWy1)./exp(1i*dWy0);
    dWx = angle(dWxss);
    dWy = angle(dWyss);
    dWx=PIE.utils.UnwrapPhaseBySortingReliabilityWithMask(dWx,255*ones(N));
    dWy=PIE.utils.UnwrapPhaseBySortingReliabilityWithMask(dWy,255*ones(N));
    dZ = dWx/2/pi;
%     dZ =PIE.utils.Retrieve_LP_iteration(dWx/2/pi,dWy/2/pi, sp, sp,ones(N))/scale;
%     dZs = PIE.utils.DelTilt(dZ);
%     
    %% plot
    xy = linspace(-L_nm/2000,L_nm/2000,N);
    figure(5),imagesc(xy,xy,dWy/2/pi);axis tight equal
    xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',14);title('Shearing phase');
    figure(2),imagesc(xy,xy,dZ);colorbar;axis tight equal
    xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',14);title('Reconstructed phase');drawnow;
    Ex_phaNs = Ex_phaN;
    % Ex_phaNs(Ex_phaNs<-0.15)=NaN;
    Ex_phaNs0 = Ex_phaNs-mean(Ex_phaNs(~isnan(Ex_phaNs)));
    nShift = round(lambda_um/T_um*Lo_um*1000/(L_nm/N));
    Ex_pha_shift = circshift(Ex_phaNs0,[0,nShift])-circshift(Ex_phaNs0,[0,-nShift]);
%     figure(3),imagesc(xy,xy,Ex_pha_shift);colorbar;axis tight equal
%     xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',14);title('Original phase');
    residual = dZ-Ex_phaN;
    residual =residual -mean(residual(:));
    res_crop = residual(30:230,30:230);
    
    RMS(j) = std(res_crop(:));
    mask_s1 = zeros(N);
    mask_s2 = zeros(N);
    rs = 15;

    shift = 26;
    mask_s1(N/2+1-rs:N/2+rs,N/2+1-rs-nShift:N/2+rs-nShift)=1;
    mask_s2(N/2+1-rs:N/2+rs,N/2+1-rs+nShift:N/2+rs+nShift)=1;
    dWs_PS(j) = abs((mean(dZ(mask_s1==1))-mean(dZ(mask_s2==1)))/2-0.346);
    
    %    residual(abs(residual)>13*std(residual(:)))=0;
%     figure(4),imagesc(xy,xy,residual);colorbar;axis tight equal
%     xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',14);title('Residual error');
end
figure(6),h=plot(aberAmp,dWs_PS);
xlabel('RMS of input aberations / waves');ylabel('Phase error / waves');set(gca,'fontSize',14);
set(h,'lineWidth',2);

