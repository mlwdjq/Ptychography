% Extracts the phase of the X and Y difference wavefronts from a single
% interferomgram using Fourier domain analysis and unwraps phase
% 
% @param {double} dImage                - interferogram
% @param {uint8} u8FilterType           - Filter type, 1: Gauss, 2: Hanning
% @param {double} dFilterWidth          - Filter characteristic width as a
%                                         ratio of f0
% @param {double} dMask                 - Wavefront domain
% 
% @return {double} dWx                  - X difference phase
% @return {double} dWy                  - Y difference phase

function [dWx, dWy, dWxNoTilt, dWyNoTilt, dWxUnwrapped, dWyUnwrapped,z1,CCDrot] = extractPhaseFD(dImage, u8FilterType, ...
                        dFilterWidth, u8UnwrapAlgorithm, dMask,dMask2, NA,T_um,z2_mm,detSize)
[sr,sc]=size(dImage);
absfimgx=zeros(sr,sc);
absfimgy=zeros(sr,sc);
% if isempty(dMask)
%     [x,y]=meshgrid(linspace(-1,1,256));
%     dMask=pinhole(256);
%     dMask(x.^2+y.^2>0.8^2)=0;
% end
%phase extraction
fimg = fftshift(fft2(dImage));
l=15;
absfimgx(:,round(sc/2)+l:sc)=abs(fimg(:,round(sc/2)+l:sc));
absfimgy(round(sr/2)+l:sr,:)=abs(fimg(round(sr/2)+l:sr,:));
[Xcenter(1),Xcenter(2)]=find(absfimgx==max(absfimgx(:)));
[Ycenter(1),Ycenter(2)]=find(absfimgy==max(absfimgy(:)));
f0=sqrt((Xcenter(1)-sr/2).^2+(Xcenter(2)-sc/2).^2);

%FilterType
switch u8FilterType
    case 1
        Xwin=gaussfilt2d(ones(sr,sc), Xcenter, dFilterWidth*f0);
        Ywin=gaussfilt2d(ones(sr,sc), Ycenter, dFilterWidth*f0);
    case 2
        Xwin=hanning2d(ones(sr,sc), Xcenter, dFilterWidth*f0);
        Ywin=hanning2d(ones(sr,sc), Ycenter, dFilterWidth*f0);
end

XfilteredPeak=fimg.*Xwin;
YfilteredPeak=fimg.*Ywin;
XrealSpaceFirstOrder=ifft2(ifftshift(XfilteredPeak));
YrealSpaceFirstOrder=ifft2(ifftshift(YfilteredPeak));
dWxUnwrapped=atan2(imag(XrealSpaceFirstOrder),real(XrealSpaceFirstOrder));
dWyUnwrapped=atan2(imag(YrealSpaceFirstOrder),real(YrealSpaceFirstOrder));
   
%unwrap phase
switch u8UnwrapAlgorithm
    case 1
        dMask(dMask~=0)=255;
        if all(size(dWxUnwrapped)==size(dMask))&&all(size(dWyUnwrapped)==size(dMask))&&all(all(dMask==0|dMask==255))&&~isempty(dMask)
            dWx=lsianalyze.utils.UnwrapPhaseBySortingReliabilityWithMask(dWxUnwrapped,dMask);
            dWy=lsianalyze.utils.UnwrapPhaseBySortingReliabilityWithMask(dWyUnwrapped,dMask);
        else
            msgbox('Please check Shearing wavefront and Mask', 'Error');
            return;
        end
    case 2
        
end
dWx=dWx-mean(dWx(dMask~=0&dMask2~=0));
dWy=dWy-mean(dWy(dMask~=0&dMask2~=0));
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