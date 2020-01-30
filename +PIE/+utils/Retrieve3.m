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

function dZ = Retrieve3(dPhiX, dPhiY, dSx, dSy, dMask)

% Scale sheared wavefronts by shear maps to create proper difference
% functions
p = dPhiX./dSx;
q = dPhiY./dSy;

% Smooth boundary conditions by tiling flipped copies of sheared wavefronts
[hei, wid] = size(p);
N=hei;                                                                                                                                                                                    
p1 =  flipud(p);
p2 = -fliplr(p);
q1 = -flipud(q);
q2 = fliplr(q);

% p = [-rot90(p,2) p1;p2 p];
% q = [-rot90(q,2) q1;q2 q];
ps=zeros(2*hei+1,2*wid+1);
qs=zeros(2*hei+1,2*wid+1);
ps(1:N,1:N)=-rot90(p,2);
ps(N+2:end,1:N)=p2;
ps(1:N,N+2:end)=p1;
ps(N+2:end,N+2:end)=p;
qs(1:N,1:N)=-rot90(q,2);
qs(N+2:end,1:N)=q2;
qs(1:N,N+2:end)=q1;
qs(N+2:end,N+2:end)=q;
% set up cdft coords:

% nr = 0:N-1;
% nc = 0:N-1;
% kr = (0:N-1)/N;
% kc = (0:N-1)/N;
[kr, kc] = meshgrid((0:N-1)/N);
% Transform sheared wavefronts and divide by frequency variable
P = fftshift(fft2(ps));%.* exp(-2i*pi  * (kr * 0.5 + kc * 0.5) ));
Q = fftshift(fft2(qs));%.* exp(-2i*pi  * (kr * 0.5 + kc * 0.5) ));
[h2,w2] = size(P);

[U,V] = meshgrid(linspace(-1/2,1/2,w2),linspace(-1/2,1/2,h2));

t1 = (P+1i*Q)./(U+1i*V)./2/pi/1i;
t1(isinf(t1)) = 0;
t1(isnan(t1))=0;
% Inverse transform and mask result
t2 = ifft2(ifftshift(t1)); 
dZ = real(t2);
dZ = dZ(hei+2:2*hei+1,wid+2:2*wid+1);   
dZ = (dZ).*dMask;
dZ(dMask==0)=NaN;

