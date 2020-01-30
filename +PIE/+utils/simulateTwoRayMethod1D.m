

% dPhShifts: [N*M x 2] array of phase shifts in radians

function [wInt, Wx, Wy] = simulateTwoRayMethod1D(...
        N, lambda_um, NA, T_um, z1_um, z2_mm, alpha, gamma,...
        zernCouples, dMaxPhoton, detectorSize,CenterObstruction, SelectMask,dPhShifts, isX) 

% Set orthogonal tilt to 0

beta    = 0;
eta     = 0;

z2_um       = z2_mm * 1e3;
halfD       = detectorSize/2*1e3;%tan(asin(NA))*(z2_um);
xIdx        = linspace(-halfD, halfD, N);
yIdx        = linspace(-halfD, halfD, N);
[X, Y]      = meshgrid(xIdx,yIdx);

Z = X*sin(eta)+Y*sin(gamma);
X = X * cos(eta);
Y = Y * cos(gamma);


% Computes the opd from +1 and 0 grating orders
if z1_um>=0
    [opds2x, opds2y] = lsianalyze.utils.twoRaySystemP(X, Y, Z,z1_um, z2_um-z1_um, lambda_um, T_um, [0, 1, 0], ...
            beta, alpha, zernCouples', NA);
    [opds2sx, opds2sy] = lsianalyze.utils.twoRaySystemP(X, Y, Z, z1_um, z2_um-z1_um, lambda_um, T_um, [0, -1, 0], ...
            beta, alpha, zernCouples', NA);
else
    [opds2x, opds2y] = lsianalyze.utils.twoRaySystem(X, Y, Z,-z1_um, z2_um-z1_um, lambda_um, T_um, [0, 1, 0], ...
            beta, alpha, zernCouples', NA);
    [opds2sx, opds2sy] = lsianalyze.utils.twoRaySystem(X, Y, Z, -z1_um, z2_um-z1_um, lambda_um, T_um, [0, -1, 0], ...
            beta, alpha, zernCouples', NA);
end
Wx=(-opds2x+ opds2sx)/(2*lambda_um); % waves
Wy=(-opds2y+ opds2sy)/(2*lambda_um); % waves     
% % Computes the opd from -1 and 1 grating orders
% [opds2x11, opds2y11] = lsianalyze.utils.twoRaySystem(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, [-1, 1, 0], ...
%         alpha, beta, zernCouples', NA);
% Computes the opd from -1 and 1 grating orders    
% [opds2sx, opds2sy] = lsianalyze.utils.twoRaySystem(X, Y, -Z, z1_um, z2_um, lambda_um, T_um, [0, -1, 1], ...
%         alpha, beta, zernCouples', NA); 
%2D
% A=[0.4,0.4,0.16,0.4,0.4,0.16,0.16,0.16,0.16,0.16];
% orders={[0, 0; 1, 0],[0, 0; -1, 0],[-1, 0; 1, 0],...
%     [0, 0; 0, 1],[0, 0; 0, -1],[0, -1; 0, 1],...
%     [-1, 0; 0, -1],[-1, 0; 0, 1],[0, -1; 1, 0],[0, 1; 1, 0]};
%1D
%orders amplitudes 
a0=sqrt(0.2689);a1=sqrt(0.1383);a2=sqrt(0.0163);
a_1=a1;
a_2=a2;
A=[a0*a1,a0*a_1,a1*a2,a_1*a_2];
if isX
orders={[1, 0; 0, 0],[-1, 0; 0, 0],[2, 0; 1, 0],[-2, 0; -1, 0]};
else
orders={[0, 1; 0, 0],[0, -1; 0, 0],[0, 2; 0, 1],[0, -2; 0, -1]};
end
%wInt=(a0^2+a1^2+a_1^2+a2^2+a_2^2)/2+1;
wInt=5*sum(A);
if z1_um>=0
    for i=1:length(orders)
        w=2*pi/lambda_um*lsianalyze.utils.twoRaySystemP2D(X, Y, Z,z1_um, z2_um-z1_um, lambda_um, T_um, ...
        orders{i}, beta, alpha, zernCouples', NA,dPhShifts/2/pi*lambda_um);
        wInt=wInt+A(i)*cos(w);
    end
else
    for i=1:length(orders)
        w=2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, Z,-z1_um, z2_um-z1_um, lambda_um, T_um, ...
        orders{i}, beta, alpha, zernCouples', NA,dPhShifts/2/pi*lambda_um);
        wInt=wInt+A(i)*cos(w);
    end
end
% w0010 = 2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, ...
%     [0, 0; 1, 0], alpha, beta, zernCouples', NA,dPhShifts/2/pi*lambda_um);
% w0001 = 2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, ...
%     [0, 0; 0, 1], alpha, beta, zernCouples', NA,dPhShifts/2/pi*lambda_um);
% w00_10 = 2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, ...
%     [0, 0; -1, 0], alpha, beta, zernCouples', NA,dPhShifts/2/pi*lambda_um);
% w000_1 = 2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, ...
%     [0, 0; 0, -1], alpha, beta, zernCouples', NA,dPhShifts/2/pi*lambda_um);
% w_1010 = 2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, ...
%     [-1, 0; 1, 0], alpha, beta, zernCouples', NA,dPhShifts/2/pi*lambda_um);
% w0_101 = 2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, ...
%     [0, -1; 0, 1], alpha, beta, zernCouples', NA,dPhShifts/2/pi*lambda_um);
% w_100_1 = 2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, ...
%     [-1, 0; 0, -1], alpha, beta, zernCouples', NA,dPhShifts/2/pi*lambda_um);
% w_1001 = 2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, ...
%     [-1, 0; 0, 1], alpha, beta, zernCouples', NA,dPhShifts/2/pi*lambda_um);
% w0_110 = 2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, ...
%     [0, -1; 1, 0], alpha, beta, zernCouples', NA,dPhShifts/2/pi*lambda_um);
% w0110 = 2*pi/lambda_um*lsianalyze.utils.twoRaySystem2D(X, Y, -Z,z1_um, z2_um, lambda_um, T_um, ...
%     [0, 1; 1, 0], alpha, beta, zernCouples', NA,dPhShifts/2/pi*lambda_um);
%2D
% wInt=4+0.4*cos(w0010)+0.4*cos(w0001)+0.4*cos(w00_10)+0.4*cos(w000_1)+...
%     0.16*cos(w_1010)+0.16*cos(w0_101)+...
%     0.16*cos(w_100_1)+0.16*cos(w_1001)+0.16*cos(w0_110)+0.16*cos(w0110);
%1D
%wInt=2+0.65*cos(w0010)+0.65*cos(w00_10)+0.42*cos(w_1010);
% Generate interferogram:
% dPhX = dPhShifts(1);
% dPhY = dPhShifts(2);
% 
% %wx = exp(2i*pi/lambda_um*opds2x)*exp(1i*dPhX) + exp(2i*pi/lambda_um*opds2sx)*exp(-1i*dPhX) ;
% %wy = exp(2i*pi/lambda_um*opds2y)*exp(1i*dPhY) + exp(2i*pi/lambda_um*opds2sy)*exp(-1i*dPhY) ;
% % wx11= exp(2i*pi/lambda_um*opds2x11)*exp(2i*dPhX);
% % wy11= exp(2i*pi/lambda_um*opds2y11)*exp(2i*dPhY);
% % wx=0.4*wx+0.16*wx11;
% % wy=0.4*wy+0.16*wy11;
% wx01=2*pi/lambda_um*opds2x+dPhX;
% wx0_1=2*pi/lambda_um*opds2sx-dPhX;
% wy01=2*pi/lambda_um*opds2y+dPhY;
% wy0_1=2*pi/lambda_um*opds2sy-dPhY;
% wx_11=2*pi/lambda_um*opds2x11+2*dPhX;
% wy_11=2*pi/lambda_um*opds2y11+2*dPhY;
    

%U = 2+ 1/pi*wx + 1/pi*wy;

% Normalize to high contrast CCD
%wInt = abs(U).^2;
wInt = (wInt/5/sum(A) * dMaxPhoton);
wInt = wInt + sqrt(wInt).*randn(size(wInt));
wInt = abs(wInt);

% crop to NA:
detRat = (2*tan(asin(NA))*(z2_um)/1000/detectorSize);
ObsRat = (2*tan(asin(NA*CenterObstruction))*(z2_um)/1000/detectorSize);
mask=pinhole(round(detRat * N), N, N)-pinhole(round(ObsRat * N), N, N);
%remove lines
switch SelectMask
    case 2
        LineWidth=4;
        HalfAngle=45;
        [x,y]=meshgrid(linspace(-1,1,N));
        mask((abs(x-tan(HalfAngle/180*pi)*y)/sqrt(1+tan(HalfAngle/180*pi)^2)<LineWidth*2/N|...
            abs(x+tan(HalfAngle/180*pi)*y)/sqrt(1+tan(HalfAngle/180*pi)^2)<LineWidth*2/N)&y>=0)=0;
end
wInt = wInt.*mask;
