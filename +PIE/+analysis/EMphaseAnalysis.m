
%% EUV mask data
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

%% QWLSI mask data
data1=load('aerial-images_2020-07-24 13.07.15.mat');
data2=load('aerial-images_2020-07-25 18.03.02.mat');
Ecx = data1.Ex_amp.*exp(1i*data1.Ex_pha)./(data2.Ex_amp.*exp(1i*data2.Ex_pha));
Ecy = data1.Ey_amp.*exp(1i*data1.Ey_pha)./(data2.Ey_amp.*exp(1i*data2.Ey_pha));
phx = angle(Ecx);
phy = angle(Ecy);

for i=1:20
%    figure(3),imagesc(phx(:,:,i)),colorbar;
%     figure(4),imagesc(phy(:,:,i)),colorbar;
%    pause(1);
   amp(i)= mean(mean(data1.Ex_amp(120:137,120:137,i))).^2+mean(mean(data1.Ey_amp(120:137,120:137,i))).^2;
   ph90(i)=mean(mean(phy(120:137,120:137,i)))-mean(mean(phy(31:80,31:80,i)));
end

% figure(2),plot(ph90)
% plot phase difference
figure(3),h=plot(atand(thetaX),ph90(1:10)/2/pi,'-o');box on
set(h,'lineWidth',2),hold on;
h=plot(atand(thetaY),[ph90(11:15),fliplr(ph90(11:15))]/2/pi,'-ro');box on
set(h,'lineWidth',2),
xlabel('Scanning angles/degree'),%ylabel('Scanning angle \theta_y');
ylabel('Phase difference/waves');set(gca,'fontSize',14);
% plot amp difference
figure(4),h=plot(atand(thetaX),amp(1:10),'-o');box on
set(h,'lineWidth',2),hold on;
h=plot(atand(thetaY),amp(11:20),'-o');box on
set(h,'lineWidth',2),
xlabel('Scanning angles/degree'),%ylabel('Scanning angle \theta_y');
ylabel('Intensity');set(gca,'fontSize',14);

