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
N= 5;
scanningRange_mm = 4;
dr = scanningRange_mm / (N-1);
dPos_mm = [0,0];
r=dr;
ds= dr/1.3;
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
dPos_mm =dPos_mm;

dPos_mm(N^2+1:end,:)=[];
length(dPos_mm)
figure(2),plot(dPos_mm(:,2),dPos_mm(:,1),'.');
%% segment scanning
r1 = 0.1043;
r2 = 0.2141;
dPhi =pi/8;
phi = [0:dPhi:2*pi-dPhi]';
dPos_mm = [r1*sin(phi),r1*cos(phi);r2*sin(phi),r2*cos(phi)];
figure(2),plot(dPos_mm(:,2),dPos_mm(:,1),'.');

%% 3D circular scanning
N= 9;
dz_mm = [-0.004,0,0.004]; 
scanningRange_mm = 0.0002;
dr = scanningRange_mm / (N-1);
dPos_mm = [0,0];
r=dr;
ds= dr/1.2;
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

% dPos_mm(N^2+1:end,:)=[];
dzs_mm = ones(size(dPos_mm,1),1)*dz_mm;
dPos_mm=[[dPos_mm;dPos_mm;dPos_mm],dzs_mm(:)];

figure(2),plot3(dPos_mm(:,3),dPos_mm(:,2),dPos_mm(:,1),'.');box on

%% save probe
save('../../data/scanning/circular_3D_seg_2.mat','dPos_mm');