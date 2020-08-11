%% this script is used to analyze the incident aberration effects of QWLSI by calibrate the null test result


%% load data
load('aerial-images_2020-08-08 12.06.13-QWSLI-Null.mat');
[N,~,nInt] = size(Ex_amp);

%% define aberration
noAber = 0;
aberAmp = linspace(0,1,101);
RMS = zeros(1, length(aberAmp));
Nz = 15;
[x,y] = meshgrid(linspace(-1.5,1.5,N));
[th,r] = cart2pol(x,y);
for k = 1:Nz
    afn = zgen([], k, 'fnr');
    temp=afn(r, th);
    Basis(:,k) = temp(:);
end
aber = zeros(N,N,Nz);

for i = 1:nInt
    dZrn = randn(Nz,1);
    temp = reshape(Basis*dZrn,N,N);
    aber(:,:,i) = temp./std(temp(:));
end
aber0 = aber;
for j = 1:length(aberAmp)
    %% load data
    load('aerial-images_2020-08-08 12.06.13-QWSLI-Null.mat');
%     [N,~,nInt] = size(Ex_amp);
aberAmp(j)
    aber = aber0*aberAmp(j);
    
    % analysis
    PIE.analysis.QWLSIIncidentAberrationAnalysis;
    dWx0 = dWx;
    dWy0 = dWy;
    
    %% load null test data
    load('aerial-images_2020-07-24 09.02.45.mat');
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
    
    dZ =PIE.utils.Retrieve_LP_iteration(dWx/2/pi,dWy/2/pi, sp, sp,ones(N))/scale;
    dZs = PIE.utils.DelTilt(dZ);
    
    %% plot
    xy = linspace(-L_nm/2000,L_nm/2000,N);
    figure(5),imagesc(xy,xy,dWy/2/pi);axis tight equal
    xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',14);title('Shearing phase');
    figure(2),imagesc(xy,xy,dZs);colorbar;axis tight equal
    xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',14);title('Reconstructed phase');
    Ex_phaNs = Ex_phaN;
    % Ex_phaNs(Ex_phaNs<-0.15)=NaN;
    Ex_phaNs0 = Ex_phaNs-mean(Ex_phaNs(~isnan(Ex_phaNs)));
    figure(3),imagesc(xy,xy,Ex_phaNs);colorbar;axis tight equal
    xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',14);title('Original phase');
    residual = dZs-Ex_phaN;
    residual =residual -mean(residual(:));
    res_crop = residual(30:230,30:230);

    RMS(j) = std(res_crop(:));
    
    %    residual(abs(residual)>13*std(residual(:)))=0;
    figure(4),imagesc(xy,xy,residual);colorbar;axis tight equal
    xlabel('x/um');ylabel('y/mm');set(gca,'fontSize',14);title('Residual error');
end
figure(6),h=plot(aberAmp,RMS);
xlabel('RMS of input aberations / waves');ylabel('RMS / waves');set(gca,'fontSize',14);
set(h,'lineWidth',2);

