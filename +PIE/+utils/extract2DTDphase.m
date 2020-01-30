function [dWx, dWy,dWxNoTilt, dWyNoTilt, dWxUnwrapped, dWyUnwrapped,z1,CCDrot] = ...
    extract2DTDphase (ceImageStack, ceMeta, dPhaseSteps, ...
    u8FourierTransformType, u8UnwrapAlgorithm, dMask,dMask2,dLowPass,u8FilterType, NA, T_um,z2_mm,detSize)


% % For now, let's stack images into a square:
% dN = sqrt(length(ceImageStack));
% ceImageStack = reshape(ceImageStack, dN, dN);
%
% % These images are rastered, so we need to flip every other one:
% for k = 1:2:dN
%     ceImageStack(:, k) =  flipud(ceImageStack(:, k));
% end
%
% % Now crop off redundant image
% ceImageStack = ceImageStack(1:dN-1,1:dN-1);
%
% dN = dN - 1;
dN = length(ceImageStack);
[dSr, dSc] = size(ceImageStack{1});
% Create x grams by y projection, and y-grams by x-projection:
dXgrams=zeros(dSr, dSc, dN);
dYgrams=zeros(dSr, dSc, dN);
PSx=zeros(dN);
PSy=zeros(dN);
Hz=zeros(dN);
dXYgram=zeros(dSr, dSc, dN, dN);

% % FilterType
% if dLowPass~=0
%     switch u8FilterType
%         case 1
%             win=gaussfilt2d(ones(dSr, dSc),[round((dSr+1)/2),round((dSc+1)/2)], dLowPass/6);
%         case 2
%             if dLowPass>min(dSr,dSc)
%                 dLowPass=min(dSr,dSc)-1;
%             end
%             win=hanning2d(ones(dSr, dSc), [round((dSr+1)/2),round((dSc+1)/2)], dLowPass);
%     end
% else
%     win=ones(dSr, dSc);
% end
for k = 1:dN
    dXgram = zeros(dSr, dSc);
    dYgram = zeros(dSr, dSc);
    for m = 1:dN
        dXgram = dXgram + ceImageStack{m, k};
        dYgram = dYgram + ceImageStack{k, m};
        dXYgram(:,:,k,m)=ceImageStack{k, m};
%                 dXgram = dXgram + real(ifft2(ifftshift(fftshift(fft2(ceImageStack{m, k})).*win)));
%                 dYgram = dYgram + real(ifft2(ifftshift(fftshift(fft2(ceImageStack{k, m})).*win)));
%                 dXYgram(:,:,k,m)=real(ifft2(ifftshift(fftshift(fft2(ceImageStack{m, k})).*win)));
        try
            PSx(m,k)=str2double(ceMeta{m,k}.DMIRetX)/T_um/5/1000*2*pi;
            PSy(m,k)=-str2double(ceMeta{m,k}.DMIRetY)/T_um/5/1000*2*pi/1.0046;
            Hz(m,k)=(str2double(ceMeta{m,k}.HSZ)-str2double(ceMeta{1,1}.HSZ));
        catch
        end
    end
    
    dXgram(dMask == 0) = 0;
    dXgrams(:,:,k) = dXgram;
    dYgram(dMask == 0) = 0;
    dYgrams(:,:,k) = dYgram;
    
end
% ShiftingError=5*Hz*tan(6/180*pi)/1000/T_um*2*pi;
try
    ShiftingError=(25*Hz*sin(6/180*pi)-5*Hz*sin(1.12/180*pi))*cos(6/180*pi)-...
        (25*Hz*cos(6/180*pi)-25*Hz*cos(1.12/180*pi))*sin(6/180*pi);
    % ShiftingErrors=(25*Hz*sin(6/180*pi)-5*Hz*sin(1.12/180*pi))*sin(6/180*pi)+(25*Hz*cos(6/180*pi)-25*Hz*cos(1.12/180*pi))*cos(6/180*pi);
    ShiftingError=ShiftingError/T_um/5/1000*2*pi/1.0046;
    PSy=PSy+ShiftingError;
catch
end
% figure(2), mesh(ShiftingError)
% % Need to reverse the sign of the y phase steps:
% dYgrams = flip(dYgrams, 3);


%phase extraction
switch u8FourierTransformType
    case 1
        xft = (fft(dXgrams, [], 3));
        yft = (fft(dYgrams, [], 3));
        Xwrap =  angle(xft(:,:,2));
        Ywrap =  angle(yft(:,:,2));
    case 2
        if PSx(1,1)==0
            if isempty(dPhaseSteps)
                [PSx,PSy]=lsianalyze.utils.FindPhaseShift(dXYgram);
            else
                [PSx,PSy]=meshgrid(dPhaseSteps(:,1),dPhaseSteps(:,2));
            end
        end
        %         x=unwrap(PSx,[],2);
        %         x=x-mean(x(:));
        %         y=unwrap(PSy,[],1);
        %         y=y-mean(y(:));
        Orders=-2:2;
        %         global px py
        %         PSx=px;
        %         PSy=py;
        [Xwrap,Ywrap]=lsianalyze.utils.DFT2forPhaseExtraction(dXYgram,PSx,PSy,Orders);
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

% top-hat filter window
iterations=100000;
if dLowPass~=0
    dWx(dMask==0)=0;
    dWy(dMask==0)=0;
    win=lsianalyze.utils.topHatWin(dSr,dLowPass);
    dWxs=real(ifft2(ifftshift(fftshift(fft2(dWx)).*win)));
    dWys=real(ifft2(ifftshift(fftshift(fft2(dWy)).*win)));
    for j=1:iterations % remove ringing error
%          RMS=std([dWxs(dMask==0)-dWx(dMask==0);dWys(dMask==0)-dWy(dMask==0)]);
%          res=dWxs-dWx;
%          res(dMask~=0)=0;
%          if j<100||(j<10000&&mod(j,100)==0)||(mod(j,1000)==0)
%               figure(2),mesh(res),title([num2str(j),' iterations']);
%               saveGif(2,'unmaskedVsIterations_diff.gif',0.01);
%          end
        dWx(dMask==0)=dWxs(dMask==0);
        dWy(dMask==0)=dWys(dMask==0);
        dWxs=real(ifft2(ifftshift(fftshift(fft2(dWx)).*win)));
        dWys=real(ifft2(ifftshift(fftshift(fft2(dWy)).*win)));
    end
    dWx=dWxs;
    dWy=dWys;
end

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

