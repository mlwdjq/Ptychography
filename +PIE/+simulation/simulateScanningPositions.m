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
dPhi =pi/4;
phi = [0:dPhi:2*pi-dPhi]';
offset = 3*tan(6/180*pi);
dPos_mm = [r1*sin(phi),r1*cos(phi)+offset;r2*sin(phi),r2*cos(phi)+offset];
figure(2),plot(dPos_mm(:,2),dPos_mm(:,1),'.');

%% round segment scanning
r1 = 0.1043;
r2 = 0.2141;
dPhi =pi/4;
phi = [0:dPhi:2*pi-dPhi]';
offset = 3*tan(6/180*pi);
dPos_mm = [r1*sin(phi),r1*cos(phi)+offset;r2*sin(phi),r2*cos(phi)+offset];
figure(2),plot(dPos_mm(:,2),dPos_mm(:,1),'.');hold on
lambda_um =13.5e-3;
k0 = 2*pi/lambda_um;
kr = k0*sin(sqrt(dPos_mm(:,1).^2+dPos_mm(:,2).^2)/pie.uiez2.get());
phi = atan2(dPos_mm(:,1),dPos_mm(:,2));
kxy = [kr.*sin(phi),kr.*cos(phi)];
dkxy = 2*pi/pie.dc_um/pie.uieRes.get();
dPosShifts = kxy./dkxy;
kxy = round(dPosShifts).*dkxy;
phis= atan2(kxy(:,1),kxy(:,2));
krs = sqrt(kxy(:,1).^2+kxy(:,2).^2);
rs = asin(krs/k0)*pie.uiez2.get();
[x,y] =pol2cart(phis,rs);
dPos_mm = [y,x];
figure(2),plot(x,y,'o');hold off

%% round segment scanning 2
r1 = 0.1043;
r2 = 0.2141;
dPhi =pi/4;
phi = [0:dPhi:2*pi-dPhi]';
offset = 3*tan(6/180*pi);
dPos_mm = [r1*sin(phi),r1*cos(phi)+offset;r2*sin(phi),r2*cos(phi)+offset];
figure(2),plot(dPos_mm(:,2),dPos_mm(:,1),'.');hold on
lambda_um =13.5e-3;
k0 = 2*pi/lambda_um;
kr = k0*sin(atan(sqrt(dPos_mm(:,1).^2+dPos_mm(:,2).^2)/pie.uieLo.get()));
phi = atan2(dPos_mm(:,1),dPos_mm(:,2));
kxy = [kr.*sin(phi),kr.*cos(phi)];
dkxy = 2*pi/pie.dc_um/pie.uieRes.get()*pie.uieMag.get();
dPosShifts = kxy./dkxy;
kxy = round(dPosShifts).*dkxy;
phis= atan2(kxy(:,1),kxy(:,2));
krs = sqrt(kxy(:,1).^2+kxy(:,2).^2);
rs = tan(asin(krs/k0))*pie.uieLo.get();
[x,y] =pol2cart(phis,rs);
dPos_mm = [y,x];
figure(2),plot(dPos_mm(:,2),dPos_mm(:,1),'o');hold off
%% 3D circular scanning
N= 9;
Nz = 11;
dz_mm = linspace(-0.004,0.004,Nz); 
scanningRange_mm = 0.0005;
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
dPos_mm=[repmat(dPos_mm,[Nz,1]),dzs_mm(:)];

figure(2),plot3(dPos_mm(:,3),dPos_mm(:,2),dPos_mm(:,1),'.');box on

%% save probe
save('../../data/scanning/CircularRoundOffsetScanningSeg16s.mat','dPos_mm');