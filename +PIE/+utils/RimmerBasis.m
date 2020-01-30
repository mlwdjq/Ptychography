function [NIx,NIy,SX,SY,reconBasis,reconBasis_downSampling]=RimmerBasis(lambda_um,NA,T_um,z1_um,z2_mm,N,alpha,gamma,Nz,detSize)
z2_um       = z2_mm * 1000;
halfD       = detSize/2*1e3;
xIdx        = linspace(-halfD, halfD, N);
yIdx        = linspace(-halfD, halfD, N);
[X, Y]      = meshgrid(xIdx,yIdx);

% Tilt angles
%alpha   = 1.12*pi/180;
beta    = 0*pi/180;
%gamma   =1.12*pi/180;
eta     = 0*pi/180;


Z = X*sin(eta)+Y*sin(gamma);
X = X*cos(eta);
Y = Y*cos(gamma);
zernCouples = [0, 0];
if z1_um>=0
    [opds2x, opds2y] = lsianalyze.utils.twoRaySystemP(X, Y, Z,z1_um, z2_um, lambda_um, T_um, [0, 1], ...
        beta, alpha, zernCouples', NA);
    [opds2sx, opds2sy] = lsianalyze.utils.twoRaySystemP(X, Y, Z, z1_um, z2_um, lambda_um, T_um, [0, -1], ...
        beta, alpha, zernCouples', NA);
else
    [opds2x, opds2y] = lsianalyze.utils.twoRaySystem(X, Y, Z,-z1_um, z2_um, lambda_um, T_um, [0, 1], ...
        beta, alpha, zernCouples', NA);
    [opds2sx, opds2sy] = lsianalyze.utils.twoRaySystem(X, Y, Z, -z1_um, z2_um, lambda_um, T_um, [0, -1], ...
        beta, alpha, zernCouples', NA);
end
opds2x=(-opds2x+ opds2sx)/2 ;
opds2y=(-opds2y+ opds2sy)/2 ;

% Null interferogram
NIx = opds2x/lambda_um;
NIy = opds2y/lambda_um;
NIx(:,:,2)=Y/halfD;
NIy(:,:,2)=-X/halfD;

[SX,SY]=lsianalyze.utils.FindShearCoordinatesa(lambda_um,z1_um,z2_um-z1_um,T_um,halfD);
SXs=SX(1:2:end, 1:2:end);
SYs=SY(1:2:end, 1:2:end);
NS=length(SX);
NSs=length(SXs);
reconBasis= zeros(NS,NS, Nz+1);
reconBasis_downSampling= zeros(NSs,NSs, Nz+1);
halfDs=tan(asin(NA))*(z2_um);
[th,r]=cart2pol(SX/halfD,SY/halfD); 
[ths,rs]=cart2pol(SXs/halfDs,SYs/halfDs); 
% Generate bases
for k = 0:Nz
    fprintf('Computing basis %d of %d\n', k, Nz);
    afn = zgen([], k, 'fnr');
    reconBasis(:,:,k+1)=afn(r, th);
    reconBasis_downSampling(:,:,k+1)=afn(rs, ths);
end
        

