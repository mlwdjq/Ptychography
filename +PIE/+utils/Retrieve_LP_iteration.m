% Reconstructs a wavefront from its sheared wavefronts using the Fourier
% Transform method. Supports variable shear
%
% @param {double} dPhiX - Sheared X wavefront
% @param {double} dPhiY - Sheared Y wavefront
% @param {double} dSx   - X shear map
% @param {double} dSy   - Y shear map
% @param {double} dMask - Wavefront mask
%
% @return {double} dZ   - Reconstructed wavefront

function dZ = Retrieve_LP_iteration(dPhiX, dPhiY, dSx, dSy,dMask,iteration)

% Scale sheared wavefronts by shear maps to create proper difference
% functions
if nargin~=6
    iteration=100;
end

[sr,sc]=size(dPhiX);
p = dPhiX;
q = dPhiY;

% Smooth boundary conditions by tiling flipped copies of sheared wavefronts
[hei, wid] = size(p);

p1 =  flipud(p);
p2 = -fliplr(p);
q1 = -flipud(q);
q2 = fliplr(q);

p = [-rot90(p,2) p1;p2 p];
q = [-rot90(q,2) q1;q2 q];

% Transform sheared wavefronts and divide by frequency variable
P = fftshift(fft2(p));
Q = fftshift(fft2(q));
[h2,w2] = size(P);
[U,V] = meshgrid(linspace(-1/2,1/2-1/w2,w2),linspace(-1/2,1/2-1/h2,h2));

% generate frequency response 
dx=-1i*sin(2*pi*dSx*(sc-1).*U);
dy=-1i*sin(2*pi*dSy*(sr-1).*V);
dxy=dx+1i*dy;

% apply low-pass filter
cutoff=0.95;
LP=double(abs(dSx*(sc-1).*U)<0.5*cutoff&abs(dSy*(sr-1).*V)<0.5*cutoff);
%  LP=lsianalyze.utils.topHatCosWin(2*sr,1.9/dSx,2/dSx);
% LP=LP&LP2;
% LP=double(~(U.^2+V.^2>0.1^2&abs(dxy)<0.1));
t1 = (P+1i*Q)./dxy;
t1(isinf(t1)) = 0;
t1(isnan(t1))=0;
% remove spectrums outside the band limit
t1=t1.*LP;
t2 = ifft2(ifftshift(t1));
dZ = real(t2);

% do iteration to remove edge effect
dMasks=[dMask,dMask;dMask,dMask];
for j=1:iteration
    dWx=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dx)));
    dWy=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dy)));
    p(dMasks==0)=dWx(dMasks==0);
    q(dMasks==0)=dWy(dMasks==0);
    P = fftshift(fft2(p));
    Q = fftshift(fft2(q));
    t1 = (P+1i*Q)./dxy;
    t1(isinf(t1)) = 0;
    t1(isnan(t1))=0;
    t1=t1.*LP;
    % Inverse transform and mask result
    t2 = ifft2(ifftshift(t1));
    dZ = real(t2);
end
dZ = dZ(hei+1:2*hei,wid+1:2*wid);
dZ = (dZ).*dMask;
dZ(dMask==0)=NaN;