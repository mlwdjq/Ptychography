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

function dZ = Retrieve_LP_iteration3(dPhiX, dPhiY, dSx, dSy,dMask,iteration)

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

factor=3/2;
ps=p*factor;
qs=p*factor;
% Transform sheared wavefronts and divide by frequency variable
P = fftshift(fft2(p));
Q = fftshift(fft2(q));
Ps = fftshift(fft2(ps));
Qs = fftshift(fft2(qs));
[h2,w2] = size(P);

[U,V] = meshgrid(linspace(-1/2,1/2-1/w2,w2),linspace(-1/2,1/2-1/h2,h2));

dx=-1i*sin(2*pi*dSx*(sc-1).*U);
dy=-1i*sin(2*pi*dSy*(sr-1).*V);
dxs=-1i*sin(2*pi*dSx*factor*(sc-1).*U);
dys=-1i*sin(2*pi*dSy*factor*(sr-1).*V);
dxy=(dx+1i*dy);
dxys=(dxs+1i*dys);
t1 = (P+1i*Q)./dxy;
t1s = (Ps+1i*Qs)./dxys;
flag=abs(dxy)<0.1&(U.^2+V.^2)>0.2^2;
t1(flag)=t1s(flag);
t1(isinf(t1)) = 0;
t1(isnan(t1))=0;
t2 = ifft2(ifftshift(t1));
dZ = real(t2);

dMasks=[dMask,dMask;dMask,dMask];

for j=1:iteration
    dWx=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dx)));
    dWy=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dy)));
    dWxs=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dxs)));
    dWys=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dys)));
    p(dMasks==0)=dWx(dMasks==0);
    q(dMasks==0)=dWy(dMasks==0);
    ps=dWxs;
    qs=dWys;
    P = fftshift(fft2(p));
    Q = fftshift(fft2(q));
    Ps = fftshift(fft2(ps));
    Qs = fftshift(fft2(qs));
    t1 = (P+1i*Q)./(dxy);
    t1s = (Ps+1i*Qs)./dxys;
    t1(flag)=t1s(flag);
    t1(isinf(t1)) = 0;
    t1(isnan(t1))=0;
    t2 = ifft2(ifftshift(t1));
    dZ = real(t2);
end

dZ = dZ(hei+1:2*hei,wid+1:2*wid);
% dZ = (dZ).*dMask;
% dZ(dMask==0)=NaN;