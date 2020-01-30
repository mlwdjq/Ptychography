% Computes interferogram based on 2-ray method

function [wInt, Wx, Wy] = simulateTwoRayMethod(...
        N, lambda_um, NA, T_um, z1_um, z2_mm, alpha, gamma,...
        zernCouples, dMaxPhoton, detectorSize,CenterObstruction, SelectMask,dPhShifts) 

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


% Optionally compute the 1D phase?
if nargout == 3
    % Computes the opd from +1 and 0 grating orders
    if z1_um>=0
        [opds2x, opds2y] = lsianalyze.utils.twoRaySystemP_rot2D(X, Y, Z,z1_um, z2_um-z1_um, lambda_um, T_um, [0, 1, 0], ...
                beta, alpha, zernCouples', NA);

        [opds2sx, opds2sy] = lsianalyze.utils.twoRaySystemP_rot2D(X, Y, Z, z1_um, z2_um-z1_um, lambda_um, T_um, [0, -1, 0], ...
                beta, alpha, zernCouples', NA);

    else
        [opds2x, opds2y] = lsianalyze.utils.twoRaySystem(X, Y, Z,-z1_um, z2_um-z1_um, lambda_um, T_um, [0, 1, 0], ...
                beta, alpha, zernCouples', NA);
        [opds2sx, opds2sy] = lsianalyze.utils.twoRaySystem(X, Y, Z, -z1_um, z2_um-z1_um, lambda_um, T_um, [0, -1, 0], ...
                beta, alpha, zernCouples', NA);
    end
    Wx=(-opds2x+ opds2sx)/(2*lambda_um); % waves
    Wy=(-opds2y+ opds2sy)/(2*lambda_um); % waves   
end



% Simulate a grating by computing effective OPD of all orders
a00=sqrt(0.2353);a10=sqrt(0.0825);a11=sqrt(0.015);a01=sqrt(0.0825);a1_1=sqrt(0.015);
a_10=a10;
a_1_1=a11;
a0_1=a01;
a_11=a1_1;

A=[a00*a10,a00*a_10,a00*a01,a00*a0_1,a01*a11,a10*a11,...
    a_10*a_1_1,a0_1*a_1_1,a10*a1_1,a0_1*a1_1,a_10*a_11,a01*a_11,...
    a01*a10,a01*a_10,a0_1*a10,a0_1*a_10];

orders={[0, 0; 1, 0],[0, 0; -1, 0],[0, 0; 0, 1],...
    [0, 0; 0, -1],[0, 1; 1, 1],[1, 0; 1, 1],...
    [-1, 0; -1, -1],[0, -1; -1, -1],[1, 0; 1, -1],...
    [0, -1; 1, -1],[-1, 0; -1, 1],[0, 1; -1, 1],...
    [0, 1; 1, 0],[0, 1; -1, 0],[0, -1; 1, 0],[0, -1; -1, 0]};


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
