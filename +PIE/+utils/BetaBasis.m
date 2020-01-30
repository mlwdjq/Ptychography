function [betaBasis,NI]=BetaBasis(lambda_um,NA,T_um,z1_um,z2_mm,N,alpha,gamma,Nz,dDetectorSize)
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


for k = 0:Nz

    if k == 0
        zernCouples = [0, 0];
 
    else
        zernCouples = [k, 1];
    end

    fprintf('Computing basis %d of %d\n', k, Nz);
    % Computes the opd from +1 and 0 grating orders
    if z1_um>=0
    [opds2x, opds2y] = lsianalyze.utils.twoRaySystemP(X, Y, Z,z1_um, z2_um-z1_um, lambda_um, T_um, [0, 1], ...
            beta, alpha, zernCouples', NA);
    [opds2sx, opds2sy] = lsianalyze.utils.twoRaySystemP(X, Y, Z, z1_um, z2_um-z1_um, lambda_um, T_um, [0, -1], ...
            beta, alpha, zernCouples', NA);
    else
    [opds2x, opds2y] = lsianalyze.utils.twoRaySystem(X, Y, Z,-z1_um, z2_um-z1_um, lambda_um, T_um, [0, 1], ...
            beta, alpha, zernCouples', NA);
    [opds2sx, opds2sy] = lsianalyze.utils.twoRaySystem(X, Y, Z, -z1_um, z2_um-z1_um, lambda_um, T_um, [0, -1], ...
            beta, alpha, zernCouples', NA);
    end
        
        
    opds2x=(-opds2x+ opds2sx)/2 ;
    opds2y=(-opds2y+ opds2sy)/2 ;
%     tempx=opds2x/lambda_um;
%     tempx(mask==0)=NaN;
    
        

   % Generate bases
    if k == 0
        NI = [opds2x ; opds2y]/lambda_um;  % Null interferogram
    else
        betaBasis(:,:,k + 1) =  [opds2x; opds2y]/lambda_um - NI;
    end
    
    
        
end

