function [betaBasis,NI,reconBasis]=FourierBasis_deferredZ(lambda_um,NA,T_um,z1_um,z2_mm,N,alpha,gamma,Nz,dDetectorSize)
z2_um       = z2_mm * 1000;
%halfD       = tan(asin(NA))*(z2_um);
halfD=dDetectorSize/2*1e3;
NAs=NA*halfD/(tan(asin(NA))*z2_um);
xIdxs        = z2_um*tan(asin(NAs*linspace(-1, 1, N)));
yIdxs        = z2_um*tan(asin(NAs*linspace(-1, 1, N)));
[X, Y]      = meshgrid(xIdxs,yIdxs);

% Tilt angles
%alpha   = 1.12*pi/180;
beta    = 0*pi/180;
%gamma   =1.12*pi/180;
eta     = 0*pi/180;


Z = X*sin(eta)+Y*sin(gamma);
X = X*cos(eta);
Y = Y*cos(gamma);

% set aberrations:
%Nz = 37; % number of zernikes
betaBasis = zeros(2*N,N, Nz + 1); % add one vector for rotation

% Now generate a rotation vector

betaBasis(:,:,1) = [Y/halfD; -X/halfD];

% Zernike basis for reconstructed wavefront
% xn=sin(atan(X./(z2_um+Z)))/NA;
% yn=sin(atan(Y./(z2_um+Z)))/NA;
% 
% [th,r]=cart2pol(xn,yn); % transfer to pupil coordinates 
th=atan2(Y,X);
r=sin(atan(sqrt(X.^2+Y.^2)./(z2_um+Z)))/NA;

reconBasis= zeros(N,N, Nz+1);



% Computes the opd from +1 and 0 grating orders using function method
[opdx1, pprx1, ppthx1, pmrx1, pmthx1] = lsianalyze.utils.twoRaySystemP2D_deferredZ_coupled(X, Y, Z,z1_um, z2_um-z1_um, lambda_um, T_um, ...
    [0,0;1,0], beta, alpha, [], NA,0);
[opdx2, pprx2, ppthx2, pmrx2, pmthx2] = lsianalyze.utils.twoRaySystemP2D_deferredZ_coupled(X, Y, Z,z1_um, z2_um-z1_um, lambda_um, T_um, ...
    [0,0;-1,0], beta, alpha, [], NA,0);
[opdy1, ppry1, ppthy1, pmry1, pmthy1] = lsianalyze.utils.twoRaySystemP2D_deferredZ_coupled(X, Y, Z,z1_um, z2_um-z1_um, lambda_um, T_um, ...
    [0,0;0,1], beta, alpha, [], NA,0);
[opdy2, ppry2, ppthy2, pmry2, pmthy2] = lsianalyze.utils.twoRaySystemP2D_deferredZ_coupled(X, Y, Z,z1_um, z2_um-z1_um, lambda_um, T_um, ...
    [0,0;0,-1], beta, alpha, [], NA,0);

opds2x=(-opdx1+ opdx2)/2 ;
opds2y=(-opdy1+ opdy2)/2 ;

NI = [opds2x ; opds2y]/lambda_um;  % Null interferogram

% Generate bases
for k = 1:Nz
    fprintf('Computing basis %d of %d\n', k, Nz);
    zfn = zgen([], k, 'fnr');
    bx1=zfn(pprx1, ppthx1) - zfn(pmrx1, pmthx1);
    bx2=zfn(pprx2, ppthx2) - zfn(pmrx2, pmthx2);
    by1=zfn(ppry1, ppthy1) - zfn(pmry1, pmthy1);
    by2=zfn(ppry2, ppthy2) - zfn(pmry2, pmthy2);
    betaBasis(:,:,k+1) =[(-bx1+bx2)/2;(-by1+by2)/2];
    
    afn = zgen([], k, 'fnr');
    reconBasis(:,:,k)=afn(r, th);
end
afn = zgen([], 0, 'fnr');
reconBasis(:,:,Nz+1)=afn(r, th);

