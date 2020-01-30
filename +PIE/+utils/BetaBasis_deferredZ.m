function [betaBasis,NI]=BetaBasis_deferredZ(lambda_um,NA,T_um,z1_um,z2_mm,N,alpha,gamma,Nz,dDetectorSize)
z2_um       = z2_mm * 1000;
halfD       = dDetectorSize/2*1e3;
xIdx        = linspace(-halfD, halfD, N);
yIdx        = linspace(-halfD, halfD, N);
[X, Y]      = meshgrid(xIdx,yIdx);

% Tilt angles

% beta and eta???
beta    = 0*pi/180;
eta     = 0*pi/180;


Z = X*sin(eta)+Y*sin(gamma);
X = X*cos(eta);
Y = Y*cos(gamma);

% set aberrations:
%Nz = 37; % number of zernikes
betaBasis = zeros(2*N,N, Nz + 1); % add one vector for rotation

% Now generate a rotation vector

betaBasis(:,:,1) = [Y/halfD; -X/halfD];




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
end

