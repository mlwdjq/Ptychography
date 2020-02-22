%% this script is used to simulate and save scanning position

%% square scanning
N= 10;
scanningRange_mm = 0.0002;
[dm,dn] = meshgrid(0: scanningRange_mm/(N-1):scanningRange_mm);
dPos_mm  = [dm(:),dn(:)];

%% spiral scannning
N= 10;
scanningRange_mm = 0.0002;
dS = scanningRange_mm / (N-1);
r0 = 0;
theta0 = 0;
theta = linspace(0,10*pi,N^2);
a = dS/2/pi;
r = a.*theta;
x = r.*cos(theta);
y = r.*sin(theta);
dPos_mm = [y,x];
figure(2),plot(x,y,'.');

%% circular scanning
N= 101;
scanningRange_mm = 0.004;
dr = scanningRange_mm / (N-1);
dPos_mm = [0,0];
r=dr;
ds= dr/1.02;
while r<scanningRange_mm/sqrt(2)
    ns = round(2*pi*r/ds);
    theta = 0:2*pi/ns:2*pi*(1-1/ns);
    x = r.*cos(theta);
    y = r.*sin(theta);
    flag = (abs(x)>scanningRange_mm/2|abs(y)>scanningRange_mm/2);
    x(flag) = [];
    y(flag) = [];
    dPos_mm = [dPos_mm;[y;x]'];
    r = r+dr;
end
dPos_mm =dPos_mm+scanningRange_mm/2;

dPos_mm(N^2+1:end,:)=[];
length(dPos_mm)
figure(2),plot(dPos_mm(:,2),dPos_mm(:,1),'.');
%% save probe
save('../../data/scanning/circular_101.mat','dPos_mm');