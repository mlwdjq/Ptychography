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

function dZ = Retrieves(dPhiX, dPhiY, dSx, dSy, dMask)

% Scale sheared wavefronts by shear maps to create proper difference
% functions
[sr,sc]=size(dPhiX);
p = dPhiX/sr;
q = dPhiY/sc;

% Smooth boundary conditions by tiling flipped copies of sheared wavefronts
[hei, wid] = size(p);


% Transform sheared wavefronts and divide by frequency variable
P = fftshift(fft2(p));
Q = fftshift(fft2(q));
[h2,w2] = size(P);

[U,V] = meshgrid(linspace(-1/2,1/2-1/w2,w2),linspace(-1/2,1/2-1/h2,h2));

 dx=(exp(2*1i*pi*dSx.*U)-exp(-2*1i*pi*dSx.*U))/2;
 dy=(exp(2*1i*pi*dSy.*V)-exp(-2*1i*pi*dSy.*V))/2;
t1 = (P+1i*Q)./(dx+1i*dy);

% t1 = (P+1i*Q)./(U+1i*V)./2/pi/1i;
t1(isinf(t1)) = 0;
t1(isnan(t1))=0;
% Inverse transform and mask result
t2 = ifft2(ifftshift(t1)); 
dZ = real(t2);
dZ = (dZ).*dMask;
dZ(dMask==0)=NaN;

