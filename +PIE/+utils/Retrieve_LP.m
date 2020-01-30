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

function dZ = Retrieve_LP(dPhiX, dPhiY, dSx, dSy)

% Scale sheared wavefronts by shear maps to create proper difference
% functions
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

 dx=1i*sin(2*pi*dSx*(sc-1).*U);%(exp(2*1i*pi*dSx*(sc-1).*U)-exp(-2*1i*pi*dSx*(sc-1).*U))/2;
 dy=1i*sin(2*pi*dSy*(sr-1).*V);%(exp(2*1i*pi*dSy*(sr-1).*V)-exp(-2*1i*pi*dSy*(sr-1).*V))/2;
t1 = (P+1i*Q)./(dx+1i*dy);

% remove spectrums outside the band limit
t1(isinf(t1)) = 0;
t1(isnan(t1))=0;
t1(abs(dSx*(sc-1).*U)>0.5&abs(dSy*(sr-1).*V)>0.5)=0;
t1(abs(dx)<1e-6&abs(dy)<1e-6)=0;

% Inverse transform and mask result
t2 = ifft2(ifftshift(t1)); 
dZ = real(t2);

dZ = dZ(hei+1:2*hei,wid+1:2*wid);   

