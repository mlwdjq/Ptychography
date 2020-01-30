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

function [dZrn, ceBasisVectors,RMSFit] = reconstructByDerivativeBasis(dWx, dWy, ...
                        dZernOrder, dMask,dMask0,lambda_um, NA, T_um, z1_um, z2_mm,...
                        N, alpha, gamma, dDetectorSize, ceBasisNIStruct,ScaledNW)    
if ~isempty(ceBasisNIStruct)
    betaBasis = ceBasisNIStruct.ceBasisVectors.betaBasis;
    NI = ceBasisNIStruct.ceBasisVectors.NI;
    NIz1_um= ceBasisNIStruct.ceBasisVectors.z1_um;
else
    [betaBasis,NI]=lsianalyze.utils.BetaBasis_deferredZ(lambda_um,NA,T_um,z1_um,z2_mm,N,alpha,gamma,dZernOrder,dDetectorSize);
    NIz1_um=z1_um;
end

%produce betaBases
% save('betaBasesdata.mat','betaBases','NI');
% if isempty(dMask)
%     [x,y]=meshgrid(linspace(-1,1,256));
%     dMask=pinhole(256);
%     dMask(x.^2+y.^2>0.8^2)=0;
% end
ceBasisVectors=cell(1,dZernOrder+1);
dMask(dMask0==0)=0;
dMask2=[dMask; dMask];
BasisVectors=ones(sum(dMask2(:)),dZernOrder+1);
for i=1:dZernOrder + 1
    ceBasisVectors{1,i}=betaBasis(:,:,i);
    BasisVectors(:,i)=ceBasisVectors{1,i}(dMask2==1);
end
NIs=NI(dMask2==1);
% NIs=[NIs,BasisVectors(:,1)];
dW=[dWx ; dWy]/2/pi; 
dWs=dW(dMask2==1);
if ScaledNW==1
%fitting best NI and remove it 
    CoefNI=pinv(NIs)*dWs;
    fprintf('NI coeff: %0.3f, with rel z1 factor: %0.3f\n', CoefNI, z1_um./NIz1_um);
else
    % Try matching z3:
    dNIDecomp   = pinv(BasisVectors)*NIs;
    dNiDef      = dNIDecomp(4);
    dWvDecomp   = pinv(BasisVectors)*dWs;
    dWvDef      = dWvDecomp(4);
    CoefNI      = dWvDef / dNiDef;
    fprintf('Matching NI coeff to def: %0.3f, ', CoefNI);
end
dWs=dWs-NIs*CoefNI;

%fitting Zernike coefficients of reconstructed wavefront
dZrn = pinv(BasisVectors)*dWs;

% remove rotation:
fprintf('Removing %0.3f waves of rotation vector\n', dZrn(1));




fittedwave=BasisVectors*dZrn;
dZrn = dZrn(2:end);
% dZrn(15)=0;
% dZrn(16)=0;
% dZrn(24)=0;
%dZrn
%%
% FW=dW;
% FW(dMask2==1)=fittedwave;
% FWX=FW(1:length(FW)/2,:);
% FWY=FW(length(FW)/2+1:end,:);
% Ori=dW;
% Ori(dMask2==1)=dWs;
% OriX=Ori(1:length(Ori)/2,:);
% OriY=Ori(length(Ori)/2+1:end,:);
% figure,mesh(OriX)
% figure,mesh(OriY)
% figure,mesh(FWX)
% figure,mesh(FWY)
%%
resWs=dWs-fittedwave;
% RMSFit=std(resWs);
ceBasisVectors = struct;
ceBasisVectors.NI = NI;
ceBasisVectors.betaBasis = betaBasis;
ceBasisVectors.z1_um = z1_um;
%%
res=dMask2;
res(dMask2==1)=resWs;
res(dMask2==0)=NaN;
resdWx=res(1:N,1:N);
resdWy=res(N+1:2*N,1:N);
% filter=fspecial('average',[5 5]); %% 先定义一个滤波器
% resdWx=imfilter(resdWx,filter,'replicate'); 
% resdWy=imfilter(resdWy,filter,'replicate'); 
RMSFit=std([resdWx(~isnan(resdWx));resdWy(~isnan(resdWy))]);
% std(resdWx(~isnan(resdWx)))
% std(resdWy(~isnan(resdWy)))
% figure(5),imagesc(resdWx);
% figure(6),imagesc(resdWy);

%% fitting derivatives independently
% for i=1:dZernOrder + 1
%     temp=betaBasis(:,:,i);
%     tempx=temp(1:N,1:N);
%     tempy=temp(N+1:2*N,1:N);
%     BasisVectorsX(:,i)=tempx(dMask==1);
%     BasisVectorsY(:,i)=tempy(dMask==1);
% end
%  NIx=NI(1:N,1:N);
%  NIy=NI(N+1:2*N,1:N);
%  NIxs=NIx(dMask==1);
%   NIys=NIy(dMask==1);
%   dWxs=dWx(dMask==1);
%   dWys=dWy(dMask==1);
% CoefNIx=pinv(NIxs)*(dWxs/2/pi);        
% CoefNIy=pinv(NIys)*(dWys/2/pi);   
% dWxs=dWxs/2/pi-CoefNIx*NIxs;
% dWys=dWys/2/pi-CoefNIy*NIys;
% 
% dZrnx = pinv(BasisVectorsX)*dWxs;
% dZrny = pinv(BasisVectorsY)*dWys;
% 
% fittedwavex=BasisVectorsX*dZrnx;
% fittedwavey=BasisVectorsY*dZrny;
% dZrnx = dZrnx(2:end);
% dZrny = dZrny(2:end);
% figure(7),bar(dZrnx);
% figure(8),bar(dZrny);
% resWsx=dWxs-fittedwavex;
% resWsy=dWys-fittedwavey;
% RMSFitx=std(resWsx);
% RMSFity=std(resWsy);