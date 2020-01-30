% Reconstructs the wavefront using
%
% @param {double} dWx               - X-sheared wavefront
% @param {double} dWy               - Y-sheared wavefront
%
% @param {doulbe} dZernOrder        - Maximum Zernike basis order
%
% @param {double} dMask             - Sheared wavefront analysis domain
%
% @return {double} dZrn             - Zernike coefficients of reconstructed wavefront
% @return {cell} ceBasisVectors     - Cell array of basis vectors

function [dZ,dZrn, ceBasisVectors,RMSFit] = reconstructByFourier(dWx, dWy, ...
    dZernOrder, dMask,dMask0,lambda_um, NA, T_um, z1_um, z2_mm,...
    N, alpha, gamma, dDetectorSize, ceBasisNIStruct,ScaledNW,RTx,RTy,RD,fittingType,orthogonalization)

% load basis
if ~isempty(ceBasisNIStruct)
    betaBasis = ceBasisNIStruct.ceBasisVectors.betaBasis;
    NI = ceBasisNIStruct.ceBasisVectors.NI;
    NIz1_um= ceBasisNIStruct.ceBasisVectors.z1_um;
    reconBasis=ceBasisNIStruct.ceBasisVectors.reconBasis;
else
    [betaBasis,NI,reconBasis]=lsianalyze.utils.FourierBasis_deferredZ(lambda_um,NA,T_um,z1_um,z2_mm,N,alpha,gamma,dZernOrder,dDetectorSize);
    NIz1_um=z1_um;
end
tic,
z2_um       = z2_mm * 1000;

halfD=dDetectorSize/2*1e3;
xIdx        = linspace(-halfD, halfD, N);
yIdx        = linspace(-halfD, halfD, N);
NAs=NA*halfD/(tan(asin(NA))*z2_um);
dSx=lambda_um/T_um/NAs/2;%shearing ratio
dSy=lambda_um/T_um/NAs/2;
[X, Y]      = meshgrid(xIdx,yIdx);
xIdxs        = z2_um*tan(asin(NAs*linspace(-1, 1, N)));
yIdxs        = z2_um*tan(asin(NAs*linspace(-1, 1, N)));
[Xs, Ys]      = meshgrid(xIdxs,yIdxs);
eta     = 0*pi/180;

X = X*cos(eta); % CCD Space
Y = Y*cos(gamma);

Xs = Xs*cos(eta); % Pupil Space
Ys = Ys*cos(gamma);

switch 1
    case 1
        interpMethod='linear';
    case 2 % best case
        interpMethod='nearest';
    case 3
        interpMethod='cubic';
    case 4
        interpMethod='spline';
end

dMask(dMask0==0)=0;

if ScaledNW~=1
    % removing NI by fitting [X;Y]
    dMaskXY=dMask;
    dMaskXYs=[dMaskXY;dMaskXY];
    XY=[X;Y]/T_um/z2_um*2*pi;
    XYs=XY(dMaskXYs==1);
    rotXY=[Y;-X]/halfD;
    rotXYs=rotXY(dMaskXYs==1);
    dWxy=[dWx;dWy];
    dWxys=dWxy(dMaskXYs==1);
    CoefNI=pinv([XYs,rotXYs])*dWxys;
%     CoefNI(1)=0;
    dWxys=dWxys-[XYs,rotXYs]*CoefNI;
    dWxy(dMaskXYs==1)=dWxys;
    dWxy(dMaskXYs==0)=0;
    dWx=dWxy(1:N,:);
    dWy=dWxy(N+1:2*N,:);
    fprintf('NI coeff corresponds to a grating defocus of z1 = %0.3f um\n',CoefNI(1));
    fprintf('Removing %0.3f waves of rotation vector\n', CoefNI(2));
end

% coordinate transfrom to pupil space
dPhiX=interp2(X,Y,dWx,Xs,Ys,interpMethod);
dPhiY=interp2(X,Y,dWy,Xs,Ys,interpMethod);
dMask=interp2(X,Y,dMask,Xs,Ys,interpMethod);
dMask(dMask~=1)=0;

dMask2=[dMask; dMask];

BasisVectors=ones(sum(dMask2(:)),dZernOrder+1);
BasisVectors_ReconstructedWavefronts=ones(sum(dMask(:)),dZernOrder+1);
for i=1:dZernOrder + 1
    temp=betaBasis(:,:,i);
    BasisVectors(:,i)=temp(dMask2==1);
    temp=reconBasis(:,:,i);
    BasisVectors_ReconstructedWavefronts(:,i)=temp(dMask==1);
end

dW=[dPhiX ; dPhiY]/2/pi;
dWs=dW(dMask2==1);

if ScaledNW==1
    %fitting best NI and remove it
    NIs=NI(dMask2==1);
    CoefNI=pinv(NIs)*dWs;
    dWs=dWs-CoefNI*NIs;
    fprintf('NI coeff: %0.3f, with rel z1 factor: %0.3f\n', CoefNI, z1_um./NIz1_um);
end

switch fittingType
    case 1
        %% fitting the reconstructed wavefront to Zernike polynomials
        zerS=zeros(dZernOrder+1,1);
        dW(dMask2==1)=dWs;
        dW(dMask2~=1)=0;
        dPhiX=dW(1:N,:);
        dPhiY=dW(N+1:2*N,:);
        % phase reconstruction
        dZ =lsianalyze.utils.Retrieve_LP_iteration(dPhiX, dPhiY, dSx, dSy,dMask);
        dZs=dZ(dMask==1);
        if orthogonalization==1 %using orthogonal basis
            [orthogonalBasis,orthogonalBasisInDomain,M,delta]=lsianalyze.utils.generateArbitraryBasisFromZernike(dMask,sin(atan(Xs/z2_um))/NA,sin(atan(Ys/z2_um))/NA,dZernOrder);
            dZrn = pinv(orthogonalBasis)*dZs;
            dZrn=M'*(dZrn./delta);
            zerS(2:4)=dZrn(2:4).*[RTx;RTy;RD];
            % RMSfit
            fittedwave=orthogonalBasis*dZrn;
            res=dZs-fittedwave;
            RMSFit=std(res);
            dZs=dZs-orthogonalBasis*zerS;
            dZ(dMask==1)=dZs;
            dZrn = dZrn(2:end);
        else
            dZrn = pinv(BasisVectors_ReconstructedWavefronts)*dZs;
            zerS(1:3)=dZrn(1:3).*[RTx;RTy;RD];
            % RMSfit
            fittedwave=BasisVectors_ReconstructedWavefronts*dZrn;
            res=dZs-fittedwave;
            RMSFit=std(res);
            dZs=dZs-BasisVectors_ReconstructedWavefronts*zerS;
            dZ(dMask==1)=dZs;
            dZrn = dZrn(1:end-1);
        end
    case 2
        %% fitting the shearing wavefront to Zernike polynomials
        zerS=zeros(dZernOrder+1,1);
        dZrn = pinv(BasisVectors)*dWs;
        zerS(2:4)=dZrn(2:4).*[RTx;RTy;RD];
        
        %RMSfit
        fittedwave=BasisVectors*dZrn;
        resWs=dWs-fittedwave;
        res=dMask2;
        res(dMask2==1)=resWs;
        res(dMask2==0)=NaN;
        resdWx=res(1:N,1:N);
        resdWy=res(N+1:2*N,1:N);
        filter=fspecial('average',[5 5]); %% 先定义一个滤波器
        resdWx=imfilter(resdWx,filter,'replicate');
        resdWy=imfilter(resdWy,filter,'replicate');
        RMSFit=std([resdWx(~isnan(resdWx));resdWy(~isnan(resdWy))]);
        
        dWs=dWs-BasisVectors*zerS;
        % remove rotation:
        fprintf('Removing %0.3f waves of rotation vector\n', dZrn(1));
        dW(dMask2==1)=dWs;
        dW(dMask2~=1)=0;
        dPhiX=dW(1:N,:);
        dPhiY=dW(N+1:2*N,:);

        dZ =lsianalyze.utils.Retrieve_LP_iteration(dPhiX, dPhiY, dSx, dSy,dMask);
        dZrn = dZrn(2:end);
        
end
dZ=interp2(Xs,Ys,dZ,X,Y,'nearest');
fprintf('Fourier reconstruction took %0.3fs\n',round(toc));

%% save basis
ceBasisVectors = struct;
ceBasisVectors.NI = NI;
ceBasisVectors.betaBasis = betaBasis;
ceBasisVectors.z1_um = z1_um;
ceBasisVectors.reconBasis =reconBasis;


