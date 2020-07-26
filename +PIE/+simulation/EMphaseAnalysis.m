%% load data
data1=load('aerial-images_2020-07-23 14.26.54.mat');
data2=load('aerial-images_2020-07-23 15.31.23.mat');
Ecx = data1.Ex_amp.*exp(1i*data1.Ex_pha)./(data2.Ex_amp.*exp(1i*data2.Ex_pha));
Ecy = data1.Ey_amp.*exp(1i*data1.Ey_pha)./(data2.Ey_amp.*exp(1i*data2.Ey_pha));
phx = angle(Ecx);
phy = angle(Ecy);

for i=1:16
%    figure(3),imagesc(phx(:,:,i)),colorbar;
%     figure(4),imagesc(phy(:,:,i)),colorbar;
%    pause(1); 
   ph90(i)=mean(mean(phy(81:100,81:100,i)))-mean(mean(phy(31:50,31:50,i)));
end

% figure(2),plot(ph90)

figure(3),h=plot3(sin(atan(dPos_mm(:,1)/3))/0.0875,sin(atan((dPos_mm(:,2)-0.324)/3))/0.0875,ph90,'-o');box on
set(h,'lineWidth',2)
xlabel('u_x'),ylabel('u_y');zlabel('Phase/rad');set(gca,'fontSize',14);