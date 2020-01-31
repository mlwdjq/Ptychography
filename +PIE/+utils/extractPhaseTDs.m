% Extracts the phase of the fundamental frequency of a stack of images
% along the 3rd dimension.
% 
% @param {N x 2 cell} ceImageStack      - cell array of images.  X images are in
%                                           column 1, Y images are in column 2
% @param {N x 2 double} dPhaseSteps     - Phase steps corresponding to images
% @param {uint8} u8FourierTransformType - Transform type, 1: FFT, 2: DFT
% @param {double} dMask                 - Wavefront domain
% 
% @return {double} dWx                  - X difference phase
% @return {double} dWy                  - Y difference phase

function [dWx, dWy,dWxNoTilt, dWyNoTilt, dWxUnwrapped, dWyUnwrapped,z1,CCDrot] = extractPhaseTDs(ceImageStack, dPhaseSteps, ...
                        u8FourierTransformType, u8UnwrapAlgorithm, dMask,dMask2, NA, T_um,z2_mm,detSize)
                    
                    
                    
% For now, let's stack images into a square:
dN = length(ceImageStack);
%ceImageStack = reshape(ceImageStack, dN, dN);

% These images are rastered, so we need to flip every other one:
% for k = 2:2:dN
%     ceImageStack(:, k) =  flipud(ceImageStack(:, k));
% end

% Now crop off redundant image
% ceImageStack = ceImageStack(1:dN-1,1:dN-1);                    
                    
                    
                    
% [num,numps]=size(ceImageStack);% Interferograms
[sr,sc]=size(ceImageStack{1,1});% Sampling x
% if isempty(dMask)
%     [x,y]=meshgrid(linspace(-1,1,256));
%     dMask=pinhole(256);
%     dMask(x.^2+y.^2>0.8^2)=0;
% end
N=sum(dMask(:));
dWx=zeros(sr,sc);
dWy=zeros(sr,sc);
dWxUnwrapped=zeros(sr,sc);
dWyUnwrapped=zeros(sr,sc);
Xgrams=zeros(N,dN,dN);
Ygrams=zeros(N,dN,dN);
Xwrap=zeros(N,1);
Ywrap=zeros(N,1);

for ps=1:dN

%obtain interferograms
for i=1:dN
    Xgrams(:,i,ps)=ceImageStack{i,ps}(dMask(:)~=0);
    Ygrams(:,i,ps)=ceImageStack{ps,i}(dMask(:)~=0);
end
end



%phase extraction
switch u8FourierTransformType
    case 1
        for j=1:N
            Xg(:,:)=Xgrams(j,:,:);
            Yg(:,:)=Ygrams(j,:,:);
            fx=fft2(Xg);
            fy=fft2(Yg);
            [Ms,Ns]=size(fx);
            if Ns>1&&Ms>1
                Xwrap(j)=atan2(imag(fx(1,2)),real(fx(1,2)));
                Ywrap(j)=atan2(imag(fy(1,2)),real(fy(1,2)));
            else
                Xwrap(j)=-atan2(imag(fx(2)),real(fx(2)));
                Ywrap(j)=-atan2(imag(fy(2)),real(fy(2)));
            end
        end 
    case 2
        for j=1:N
%             fx=dft(Xgrams(j,:),dPhaseSteps(:,1)');
%             fy=dft(Ygrams(j,:),dPhaseSteps(:,2)');
            nr=[0:length(Xgrams(j,:))-1]/length(Xgrams(j,:));
            fx=lsianalyze.utils.mcdft2(Xgrams(j,:), 0, dPhaseSteps(:,1)'/2/pi*length(Xgrams(j,:)), 0, nr);
            fy=lsianalyze.utils.mcdft2(Ygrams(j,:), 0, dPhaseSteps(:,2)'/2/pi*length(Xgrams(j,:)), 0, nr);
            
            Xwrap(j)=atan2(imag(fx(2)),real(fx(2)));
            Ywrap(j)=atan2(imag(fy(2)),real(fy(2)));
        end         
end
%unwrap phase
switch u8UnwrapAlgorithm
    case 1
        dWxUnwrapped(dMask~=0)=Xwrap;
        dWyUnwrapped(dMask~=0)=Ywrap;
        dMask(dMask~=0)=255;
        if all(size(dWxUnwrapped)==size(dMask))&&all(size(dWyUnwrapped)==size(dMask))&&all(all(dMask==0|dMask==255))&&~isempty(dMask)
            dWx=lsianalyze.utils.UnwrapPhaseBySortingReliabilityWithMask(dWxUnwrapped,dMask);
            dWy=lsianalyze.utils.UnwrapPhaseBySortingReliabilityWithMask(dWyUnwrapped,dMask);
        else
            msgbox('Please check Shearing wavefront and Mask', 'Error');
            return;
        end
%     case 2
        
end

dWx=dWx-mean(dWx(dMask~=0));
dWy=dWy-mean(dWy(dMask~=0));
dWx(dMask==0)=0;
dWy(dMask==0)=0;
dWx(dMask2==0)=0;
dWy(dMask2==0)=0;
dWxNoTilt=dWx;
dWyNoTilt=dWy;
dWxNoTilt(dMask==0)=NaN;
dWyNoTilt(dMask==0)=NaN;
dWxNoTilt(dMask2==0)=NaN;
dWyNoTilt(dMask2==0)=NaN;
[dWxNoTilt,cx]=lsianalyze.utils.DelTilt(dWxNoTilt);
[dWyNoTilt,cy]=lsianalyze.utils.DelTilt(dWyNoTilt);
% estimate z1
z1=(cx(1)+cy(2))/2*T_um/detSize*2*z2_mm/2/pi;
CCDrot=(atan2(cy(1),cy(2))/pi*180-atan2(cx(2),cx(1))/pi*180)/2;
% dWxs(:,:,ps)=dWx;
% dWys(:,:,ps)=dWy;
% dWxNoTilts(:,:,ps)=dWxNoTilt;
% dWyNoTilts(:,:,ps)=dWyNoTilt;


%avarage operation 
% dWx=mean(dWxs,3);
% dWy=mean(dWys,3);
% dWxNoTilt=mean(dWxNoTilts,3);
% dWyNoTilt=mean(dWyNoTilts,3);
% z1=mean(z1s);
% CCDrot=mean(CCDrots);