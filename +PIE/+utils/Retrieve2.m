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

function dZs = Retrieve2(dPhiX, dPhiY, dSx, dSy, dMask)

% Scale sheared wavefronts by shear maps to create proper difference
% functions
p = dPhiX./dSx;
q = dPhiY./dSy;

% Smooth boundary conditions by tiling flipped copies of sheared wavefronts
%[hei, wid] = size(p);
p1=p(:);
q1=q(:);
p2=flipud(p1);
q2 = -flipud(q1);
ps = [0;p1;p2];
qs = [0;q1;q2];
% Transform sheared wavefronts and divide by frequency variable
P = (fft(ps));
Ps = fftshift(P);
Q = (fft2(qs));
Qs = fftshift(Q);
l = length(Ps);

%[Us,Vs] = meshgrid(linspace(-1/2,1/2-1/w2,w2),linspace(-1/2,1/2-1/h2,h2));
Us=linspace(-1/2,1/2,l)';
Vs=Us;

t1 = (Ps+1i*Qs)./(Us+1i*Vs)./2/pi/1i;
t1(isinf(t1)) = 0;
t1(isnan(t1))=0;

% Inverse transform and mask result
t2 = ifft2(ifftshift(t1)); 
% figure(5),imagesc(real(t2))
% figure(6),imagesc(imag(t2))
dZ = real(t2);
dZ = dZ(2:(l-1)/2+1);
dZs=zeros(size(p));
dZs(dMask==1) = dZ;
dZs(dMask==0)=NaN;

