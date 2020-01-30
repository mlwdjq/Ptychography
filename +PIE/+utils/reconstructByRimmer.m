% Reconstructs the wavefront using
%
% @param {double} dWx                   - X-sheared wavefront
% @param {double} dWy                   - Y-sheared wavefront
%
% @param {uint8} u8RimmerType           - Rimmer type:
%                                           1: 2X downsample
%                                           2: <to be determined>
%
% @param {double} dMask                 - Sheared wavefront analysis domain
%
% @return {double} dW                   - Reconstructed wavefront
% @return {sparse double} sdRimMatrix   - Sparse Rimmer reconstruction matrix (inv(T) * A)
% @param {double} X                     - x coordinate at detector plane
% @param {double} Y                     - x coordinate at detector plane
% @param {double} T_um                  - grating pitch
function [dW,dZrn,ceBasisVectors,RMSFit] = reconstructByRimmer(dWx, dWy, ...
    u8RimmerType, dMask,dMask0,lambda_um,T_um,NA,z1_um,z2_mm, N, NZernike,alpha,gamma,detSize,ceBasisNIStruct,ScaledNW)
% load basis
if ~isempty(ceBasisNIStruct)
    NIx = ceBasisNIStruct.ceBasisVectors.NIx;
    NIy = ceBasisNIStruct.ceBasisVectors.NIy;
    SX =ceBasisNIStruct.ceBasisVectors.SX;
    SY =ceBasisNIStruct.ceBasisVectors.SY;
    reconBasis =ceBasisNIStruct.ceBasisVectors.reconBasis;
    reconBasis_downSampling =ceBasisNIStruct.ceBasisVectors.reconBasis_downSampling;
else
    [NIx,NIy,SX,SY,reconBasis,reconBasis_downSampling]=lsianalyze.utils.RimmerBasis(lambda_um,NA,T_um,z1_um,z2_mm,N,alpha,gamma,NZernike,detSize);
end
tic,
z2_um       = z2_mm * 1e3;
halfD       = detSize/2*1e3;
xIdx        = linspace(-halfD, halfD, N);
yIdx        = linspace(-halfD, halfD, N);
[X, Y]      = meshgrid(xIdx,yIdx);

dMask(dMask0==0)=0;
for n=1:2
    tempx=NIx(:,:,n);
    tempy=NIy(:,:,n);
    NIs(:,n)=[tempx(dMask==1);tempy(dMask==1)];
end

%fitting best NI and remove it

if ScaledNW==0
    % removing NI and rotation by fitting [X;Y]
    dWxys=[dWx; dWy];
    dMaskXYs=[dMask;dMask];
    XY=[X;Y]/T_um/z2_um*2*pi;
    XYs=XY(dMaskXYs==1);
    rotXY=[Y;-X]/halfD;
    rotXYs=rotXY(dMaskXYs==1);
    dWxy=dWxys(dMaskXYs==1);
    CoefNI=pinv([XYs,rotXYs])*dWxy;
    dWx=dWx/2/pi-CoefNI(1)*X/T_um/z2_um-CoefNI(2)*Y/halfD;
    dWy=dWy/2/pi-CoefNI(1)*Y/T_um/z2_um-CoefNI(2)*(-X)/halfD;
    fprintf('NI coeff corresponds to a grating defocus of z1 = %0.3f um\n',CoefNI(1));
    fprintf('Removing %0.3f waves of rotation vector\n', CoefNI(2));
else
    % removing NI by fitting model
    dWxy=[dWx(dMask==1) ; dWy(dMask==1)];
    CoefNI=pinv(NIs)*dWxy;
%     fprintf('NI coeff corresponds to a grating defocus of z1 = %0.3f um\n',CoefNI(1));
    fprintf('Removing %0.3f waves of rotation vector\n', CoefNI(2));
    dWx=dWx/2/pi-CoefNI(1)*NIx(:,:,1)-CoefNI(2)*NIx(:,:,2);
    dWy=dWy/2/pi-CoefNI(1)*NIy(:,:,1)-CoefNI(2)*NIy(:,:,2);
end

dWx(dMask==0)=0;
dWy(dMask==0)=0;

[sr,sc]=size(dWx);
if sr~=sc
    msgbox('Row must be equal to column£¡', 'Error');
    return;
end

%eta=0;
% Z = X*sin(gamma)+Y*sin(eta);
% X = X*cos(gamma);
% Y = Y*cos(eta);

%rimmer reconstruction
switch u8RimmerType
    case 1
        %         [SX,SY]=lsianalyze.utils.FindShearCoordinatesa(lambda_um,z1_um,z2_um-z1_um,T_um,halfD);
        dWx=interp2(X,Y,dWx,SX,SY);
        dWy=interp2(X,Y,dWy,SX,SY);
        phx=dWx(1:2:end, 2:2:end);
        phy=dWy(2:2:end, 1:2:end);
        dMasks=interp2(X,Y,dMask,SX,SY);
        dMasks(dMasks~=1)=0;
        dMaskx=dMasks(1:2:end, 2:2:end);
        dMasky=dMasks(2:2:end, 1:2:end);
        dMasks=dMasks(1:2:end, 1:2:end);
        BasisVectors_ReconstructedWavefronts_downSampling=ones(sum(dMasks(:)),NZernike);
        for i=1:NZernike + 1
            temp=reconBasis_downSampling(:,:,i);
            BasisVectors_ReconstructedWavefronts_downSampling(:,i)=temp(dMasks==1);
        end
        if mod(length(dWx),2)==1
            phx(:,end+1) = 0;
            phy(end+1,:) = 0;
            dMaskx(:,end+1) = 0;
            dMasky(end+1,:) = 0;
        end
        % figure,imagesc(phx)
        if z1_um>=0
            px=-phx*2;
            py=-phy*2;
        else
            px=phx*2;
            py=phy*2;
        end
        [dW, T, A] = lsianalyze.utils.rimrecon3m_ds(px, py,dMaskx,dMasky);
        dWs=dW(dMasks==1);
        dZrn = pinv(BasisVectors_ReconstructedWavefronts_downSampling)*dWs;
        fittedwave=BasisVectors_ReconstructedWavefronts_downSampling*dZrn;
        dZrn = dZrn(2:end);
        dW(dMasks==0)=0;
    case 2
        %         [SX,SY]=lsianalyze.utils.FindShearCoordinates2a(lambda_um,z1_um,z2_um-z1_um,T_um,halfD);
        dWx=interp2(X,Y,dWx,SX,SY);
        dWy=interp2(X,Y,dWy,SX,SY);
        dWx=circshift(dWx, [0, -1]);
        dWy=circshift(dWy, [-1, 0]);
        dMasks=interp2(X,Y,dMask,SX,SY);
        dMasks(dMasks~=1)=0;
        BasisVectors_ReconstructedWavefronts=ones(sum(dMasks(:)),NZernike+1);
        for i=1:NZernike + 1
            temp=reconBasis(:,:,i);
            BasisVectors_ReconstructedWavefronts(:,i)=temp(dMasks==1);
        end
        dMaskx=circshift(dMasks, [0, -1]);
        dMasky=circshift(dMasks, [-1, 0]);
        if z1_um>=0
            px=-dWx*2;
            py=-dWy*2;
        else
            px=dWx*2;
            py=dWy*2;
        end
        [dW, T, A] = lsianalyze.utils.rimrecon3m( px, py,dMaskx,dMasky);
        %         halfDs=tan(asin(NA))*(z2_um);
        %         [c,fittedwave]=lsianalyze.utils.ZernikeFitWithUnevenXY(dW,dMasks,SX/halfDs,SY/halfDs);
        dWs=dW(dMasks==1);
        dZrn = pinv(BasisVectors_ReconstructedWavefronts)*dWs;
        fittedwave=BasisVectors_ReconstructedWavefronts*dZrn;
        dZrn = dZrn(2:end);
        dW(dMasks==0)=0;
end
fw=dWs-fittedwave;
RMSFit=std(fw);
fprintf('Rimmer reconstruction took %0.3fs\n',round(toc));
% sdRimMatrix=T\A;
ceBasisVectors = struct;
ceBasisVectors.NIx = NIx;
ceBasisVectors.NIy = NIy;
ceBasisVectors.SX = SX;
ceBasisVectors.SY = SY;
ceBasisVectors.reconBasis = reconBasis;
ceBasisVectors.reconBasis_downSampling = reconBasis_downSampling;