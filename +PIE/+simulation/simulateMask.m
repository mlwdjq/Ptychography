%% this script is used to simulate and save analysis mask

%% contact mask
N = length(pie.dObject);
mask = zeros(N);
mask(round(5*N/8)+1-round(N/4):round(5*N/8)+round(N/4),round(5*N/8)+1-round(N/4):round(5*N/8)+round(N/4))=1;
mask(round(5*N/8)+1-round(N/8)-2:round(5*N/8)+round(N/8)+2,round(5*N/8)+1-round(N/8)-2:round(5*N/8)+round(N/8)+2)=0;
mask(round(5*N/8)+1-round(N/8)+2:round(5*N/8)+round(N/8)-2,round(5*N/8)+1-round(N/8)+2:round(5*N/8)+round(N/8)-2)=2;
figure(2),imagesc(mask)
%% save probe
save([pie.cAppPath,'/../../data/mask/contact.mat'],'mask');