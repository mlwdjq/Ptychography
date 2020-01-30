%% this function apply iteration process to the FFT reconstruction method, which has better performance if there is an aperture
function dZ=iterationRetrieve2(dPhiX, dPhiY, dSx, dSy,dMask,iterations)
if nargin<6
    iterations=50;
end

[sr,sc]=size(dPhiX);
if ~all(dMask(:)~=0)
    for loop=1:iterations
        dZ = -lsianalyze.utils.Retrieve_LP(dPhiX, dPhiY, dSx, dSy);
        [U,V] = meshgrid(linspace(-1/2,1/2-1/sc,sc),linspace(-1/2,1/2-1/sr,sr));
        dx=1i*sin(2*pi*dSx*(sc-1).*U);
        dy=1i*sin(2*pi*dSy*(sr-1).*V);
        dWx=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dx)));
        dWy=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dy)));
        dPhiX(dMask==0)=dWx(dMask==0);
        dPhiY(dMask==0)=dWy(dMask==0);
        res=[dWx(dMask==0)-dPhiX(dMask==0);dWy(dMask==0)-dPhiX(dMask==0)];
        rms=std(res);
        if rms<1e-5
            break;
        end
    end
    dZ(dMask==0)=NaN;
else
    dZ = -lsianalyze.utils.Retrieve(dPhiX, dPhiY, dSx, dSy, dMask);
end
% figure, mesh(dPhiX)
% % do iteration to remove edge effect
% dMasks=[dMask,dMask;dMask,dMask];
% for j=1:1000
% dWx=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dx)));
% dWy=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dy)));
% p(dMasks==0)=dWx(dMasks==0);
% q(dMasks==0)=dWy(dMasks==0);
% P = fftshift(fft2(p));
% Q = fftshift(fft2(q));
% t1 = (P+1i*Q)./(dx+1i*dy);
% t1(isinf(t1)) = 0;
% t1(isnan(t1))=0;
% t1(abs(dSx*(sc-1).*U)>0.5&abs(dSy*(sr-1).*V)>0.5)=0;
% t1(abs(dx)<1e-6&abs(dy)<1e-6)=0;
% % Inverse transform and mask result
% t2 = ifft2(ifftshift(t1)); 
% dZ = real(t2);
% end
