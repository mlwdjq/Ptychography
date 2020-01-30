function [dWx, dWy,dWxNoTilt, dWyNoTilt, dWxUnwrapped, dWyUnwrapped,z1,CCDrot] = ...
    extract2x1Dphase (ceImageStack, dPhaseSteps, ...
                        u8FourierTransformType, u8UnwrapAlgorithm, dMask,dMask2,dLowPass,u8FilterType, NA, T_um,z2_mm,detSize)

                    
% % For now, let's stack images into a square:

dN = size(ceImageStack, 1);
[dSr, dSc] = size(ceImageStack{1});
% Create x grams by y projection, and y-grams by x-projection:
dXgrams=zeros(dSr, dSc, dN);
dYgrams=zeros(dSr, dSc, dN);

%FilterType
if dLowPass~=0
    switch u8FilterType
        case 1
            win=gaussfilt2d(ones(dSr, dSc),[round((dSr+1)/2),round((dSc+1)/2)], dLowPass/6);
        case 2
            if dLowPass>min(dSr,dSc)
                dLowPass=min(dSr,dSc)-1;
            end
            win=hanning2d(ones(dSr, dSc), [round((dSr+1)/2),round((dSc+1)/2)], dLowPass);
    end
else
    win=ones(dSr, dSc);
end

for k = 1:dN
%     dXgram = ceImageStack{k, 1};
%     dYgram = ceImageStack{k, 2};
    dXgram =  real(ifft2(ifftshift(fftshift(fft2(ceImageStack{k, 1})).*win)));
    dYgram =  real(ifft2(ifftshift(fftshift(fft2(ceImageStack{k, 2})).*win)));
    
    dXgram(dMask == 0) = 0;
    dXgrams(:,:,k) = dXgram;
    dYgram(dMask == 0) = 0;
    dYgrams(:,:,k) = dYgram;
    
end

%phase extraction
switch u8FourierTransformType
    case 1
        xft = (fft(dXgrams, [], 3));
        yft = (fft(dYgrams, [], 3));
        
        % figure out number of periods:
        [~, idx] = max(abs(squeeze(xft(round(200/650*dSr), round(200/650*dSr), 2:4))));
        N = idx + 1;
        fprintf('Using %d periods\n', N - 1);
        
        Xwrap =  angle(xft(:,:,N));
        Ywrap =  angle(yft(:,:,N));
end
%unwrap phase
switch u8UnwrapAlgorithm
    case 1
%         dWxUnwrapped(dMask~=0)=Xwrap;
%         dWyUnwrapped(dMask~=0)=Ywrap;

        dWxUnwrapped = Xwrap;
        dWyUnwrapped = Ywrap;
        dMask(dMask~=0)=255;
        if all(size(dWxUnwrapped)==size(dMask))&&all(size(dWyUnwrapped)==size(dMask))&&all(all(dMask==0|dMask==255))&&~isempty(dMask)
            dWx=lsianalyze.utils.UnwrapPhaseBySortingReliabilityWithMask(dWxUnwrapped,dMask);
            dWy=lsianalyze.utils.UnwrapPhaseBySortingReliabilityWithMask(dWyUnwrapped,dMask);
        else
            msgbox('Please check Shearing wavefront and Mask', 'Error');
            return;
        end
        
        
     case 2
         % Assume here that we have phase steps for each image
         
         
        
end

dWx=dWx-mean(dWx(dMask~=0));
dWy=dWy-mean(dWy(dMask~=0));
% % top hat cos filter window
% if dLowPass~=0
%     dcutoff=40;
%     win=lsianalyze.utils.topHatCosWin(dSr,dLowPass,dcutoff);
%     dWx=real(ifft2(ifftshift(fftshift(fft2(dWx)).*win)));
%     dWy=real(ifft2(ifftshift(fftshift(fft2(dWy)).*win)));
% end

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