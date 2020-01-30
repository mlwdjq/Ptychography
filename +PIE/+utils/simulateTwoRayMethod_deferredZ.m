% Computes interferogram based on 2-ray method

function wInt = simulateTwoRayMethod_deferredZ(...
    N, lambda_um, NA, T_um, z1_um, z2_mm, alpha, gamma,...
    zernCouples,zfn, dMaxPhoton, dShot2Shot, A_11_strength,...
    A_20_strength,detectorCurve, FlareLevel,dMSFN_Phase,airflow,...
    detectorSize,CenterObstruction, SelectMask,dPhShifts)

% Set orthogonal tilt to 0
beta    = 0;
eta     = 0;

z2_um       = z2_mm * 1e3;
halfD       = detectorSize/2*1e3;%tan(asin(NA))*(z2_um);
xIdx        = linspace(-halfD, halfD, N);
yIdx        = linspace(-halfD, halfD, N);
[X, Y]      = meshgrid(xIdx,yIdx);

% here we assume the detector is curved like a cylinder
CurveZ=detectorCurve*(X/halfD).^2;

Z = X*sin(eta)+Y*sin(gamma)+CurveZ;
X = X * cos(eta);
Y = Y * cos(gamma);

% define pupil corrdinates for MSFN
nMSFN=length(dMSFN_Phase);
NAs=NA*halfD/(tan(asin(NA))*z2_um);
[xPupil, yPupil]    = meshgrid(NAs/NA*linspace(-1, 1, nMSFN));

% define airflow
nAirflow = 512;
if nAirflow<N
    nAirflow = N;
end
dAirflow_Phase =   lsianalyze.utils.generateMSFN(airflow,nAirflow,1,10); % generate MSFN
[xPupil2, yPupil2]    = meshgrid(NAs/NA*linspace(-1, 1, nAirflow));
% Simulate a grating by computing effective OPD of all orders

% WZ stuff:
% a00=sqrt(0.2353);a10=sqrt(0.0825);a11=sqrt(0.015);a01=sqrt(0.0825);a1_1=sqrt(0.015);
% a_10=a10;
% a_1_1=a11;
% a0_1=a01;
% a_11=a1_1;

% A=[a00*a10,a00*a_10,a00*a01,a00*a0_1,a01*a11,a10*a11,...
%     a_10*a_1_1,a0_1*a_1_1,a10*a1_1,a0_1*a1_1,a_10*a_11,a01*a_11,...
%     a01*a10,a01*a_10,a0_1*a10,a0_1*a_10];
%
% orders={[0, 0; 1, 0],[0, 0; -1, 0],[0, 0; 0, 1],...
%     [0, 0; 0, -1],[0, 1; 1, 1],[1, 0; 1, 1],...
%     [-1, 0; -1, -1],[0, -1; -1, -1],[1, 0; 1, -1],...
%     [0, -1; 1, -1],[-1, 0; -1, 1],[0, 1; -1, 1],...
%     [0, 1; 1, 0],[0, 1; -1, 0],[0, -1; 1, 0],[0, -1; -1, 0]};

% RM 07/05/19  Changing these orders to a CB grating
a00     = 1;
a10     = 0.408;
a11     = A_11_strength;
a01     = 0.408;
a1_1    = A_11_strength;


a_10=a10;
a_1_1=a11;
a0_1=a01;
a_11=a1_1;

a_20 = A_20_strength;
a20 = A_20_strength;

A=[a00*a10,a00*a_10,a00*a01,a00*a0_1,a01*a11,a10*a11,...
    a_10*a_1_1,a0_1*a_1_1,a10*a1_1,a0_1*a1_1,a_10*a_11,a01*a_11,...
    a01*a10,a01*a_10,a0_1*a10,a0_1*a_10, a_10*a_20, a10*a20];
B=a00^2+a10^2+a_10^2+a0_1^2+a01^2+a_1_1^2+a11^2+a1_1^2+a_11^2+a_20^2+a20^2;%background intensity
% A = A/(sum(A));


% Mixing 1 and 2 orders
orders={[0, 0; 1, 0],[0, 0; -1, 0],[0, 0; 0, 1],...
    [0, 0; 0, -1],[0, 1; 1, 1],[1, 0; 1, 1],...
    [-1, 0; -1, -1],[0, -1; -1, -1],[1, 0; 1, -1],...
    [0, -1; 1, -1],[-1, 0; -1, 1],[0, 1; -1, 1],...
    [0, 1; 1, 0],[0, 1; -1, 0],[0, -1; 1, 0],[0, -1; -1, 0], ...
    [1, 0; 2,0], [-1, 0; -2, 0]};


background=(2*sum(A)+B)/(1-FlareLevel)-2*sum(A);
wInt=background;% total background intensity


tic
for i=1:length(orders)
    [opd, ppr, ppth, pmr, pmth] = lsianalyze.utils.twoRaySystemP2D_deferredZ_coupled(X, Y, Z,z1_um, z2_um-z1_um, lambda_um, T_um, ...
        orders{i}, beta, alpha, [], NA,dPhShifts/2/pi*lambda_um);
    %         opd = lsianalyze.utils.twoRaySystemP2D(X, Y, Z,z1_um, z2_um-z1_um, lambda_um, T_um, ...
    %             orders{i}, beta, alpha, [], NA,dPhShifts/2/pi*lambda_um);
    if nMSFN>0 %adding middle spatial frequency
        [px,py]=pol2cart(ppth, ppr);
        [mx,my]=pol2cart(pmth, pmr);
        diffMSFN=interp2(xPupil, yPupil,dMSFN_Phase,px,py,'nearest') - interp2(xPupil, yPupil,dMSFN_Phase,mx,my,'nearest');
        diffMSFN(isnan(diffMSFN))=0;
    else
        diffMSFN=0;
    end
    if airflow>0 %adding airflow
        [px,py]=pol2cart(ppth, ppr);
        [mx,my]=pol2cart(pmth, pmr);
        diffAirflow=interp2(xPupil2, yPupil2,dAirflow_Phase,px,py,'nearest') - interp2(xPupil2, yPupil2,dAirflow_Phase,mx,my,'nearest');
        diffAirflow(isnan(diffAirflow))=0;
    else
        diffAirflow=0;
    end
    % Add contribution of zernikes all at once:
    %     tic
    zrnOPD = lambda_um*(zfn(ppr, ppth) - zfn(pmr, pmth)+diffMSFN + diffAirflow);
    %     fprintf('Evaluating zernike lambda took %s\n', s2f(toc));
    % mesh(-opd/lambda_um-X/T_um/z2_um*z1_um)
    w=2*pi/lambda_um * (opd + zrnOPD);
    wInt=wInt+2*A(i)*cos(w);
end
wInt=wInt/(2*sum(A)+background);% normalizing intensity

fprintf('Computing order stack for single phase step took %s\n', s2f(toc));


if z1_um < 0
    
    error('This version of simulation does not support negative z1');
end


% wInt = abs(wInt);

% crop to NA:
detRat = (2*tan(asin(NA))*(z2_um)/1000/detectorSize);
ObsRat = (2*tan(asin(NA*CenterObstruction))*(z2_um)/1000/detectorSize);
if mod(N,2)==0
    mask=pinhole(2*round(detRat/2 * N), N, N)-pinhole(2*round(ObsRat/2 * N), N, N);
else
    mask=pinhole(2*ceil(detRat/2 * N)-1, N, N)-pinhole(2*floor(ObsRat/2 * N)+1, N, N);
end
%remove lines
switch SelectMask
    case 2
        LineWidth=4;
        HalfAngle=45;
        [x,y]=meshgrid(linspace(-1,1,N));
        mask((abs(x-tan(HalfAngle/180*pi)*y)/sqrt(1+tan(HalfAngle/180*pi)^2)<LineWidth*2/N|...
            abs(x+tan(HalfAngle/180*pi)*y)/sqrt(1+tan(HalfAngle/180*pi)^2)<LineWidth*2/N)&y>=0)=0;
    case 3
        if mod(N,2)==0
            mask=lsianalyze.utils.elipticalHole(2*round(detRat/2 * N),2*round(detRat * N*0.5714/2), N, N)-...
                lsianalyze.utils.elipticalHole(2*round(ObsRat/2 * N),2*round(ObsRat*0.5714/2 * N),N,N);
        else
            mask=lsianalyze.utils.elipticalHole(2*ceil(detRat/2 * N)-1,2*ceil(detRat * N*0.5714/2)-1, N, N)-...
                lsianalyze.utils.elipticalHole(2*floor(ObsRat/2 * N)+1,2*floor(ObsRat*0.5714/2 * N)+1,N,N);
        end
end

% photon noise
if (dMaxPhoton > 0)
    wInt = (wInt * dMaxPhoton);
    wInt_min = min(wInt(mask==1));
    if wInt_min<0
        wInt = wInt + sqrt(wInt-wInt_min).*randn(size(wInt));
    else
        wInt = wInt + sqrt(wInt).*randn(size(wInt));
    end
end
% Shot to shot
wInt = wInt * (1 + dShot2Shot/100 * randn(1));

wInt = wInt.*mask;
