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

function dZ = Retrieve2Dflip(dPhiX,dPhiY, dSx,dSy,  dMask)

% Scale sheared wavefronts by shear maps to create proper difference
% functions
[sr,sc]=size(dPhiX);
p = dPhiX;
q = dPhiY;
p1 =  flipud(p);
p2 = -fliplr(p);
q1 = -flipud(q);
q2 = fliplr(q);

p = [-rot90(p,2) p1;p2 p];
q = [-rot90(q,2) q1;q2 q];


% Smooth boundary conditions by tiling flipped copies of sheared wavefronts
                                                                                                                                                                                           
% p1 =  flipud(p);
% p2 = -fliplr(p);
% q1 = -flipud(q);
% q2 = fliplr(q);
% 
% p = [-rot90(p,2) p1;p2 p];
% q = [-rot90(q,2) q1;q2 q];

% Transform sheared wavefronts and divide by frequency variable
P = (fft2(p));
Q = (fft2(q));
% Q = fftshift(fft2(q));
[h2,w2] = size(P);

% U = linspace(0,1-1/w2,w2);
[U,V] = meshgrid(linspace(0,1-1/w2,w2),linspace(0,1-1/h2,h2));

 dx=2*1i*sin(2*pi*dSx.*U);
 dy=2*1i*sin(2*pi*dSy.*V);
t1 = (P+1i*Q)./(dx+1i*dy);
[kx,ky]=find(abs(dx+1i*dy)<1e-10);
for p=1:length(kx)
    for q=1:length(ky)
        if kx(p)==1|| kx(p)==2*sr||ky(q)==1|| ky(q)==2*sc
            t1(kx(p),ky(q))=0;
        else
            t1(kx(p),ky(q))=(t1(kx(p)+1,ky(q))+t1(kx(p)-1,ky(q))+t1(kx(p),ky(q)+1)+t1(kx(p),ky(q)-1))/4;
        end
    end
end
% t1 = (P+1i*Q)./(U+1i*V)./2/pi/1i;
% t1(isinf(t1)) = 0;
% t1(isnan(t1))=0;
% t1(abs(real(t1))>1e6)=0;
% t1(abs(imag(t1))>1e6)=0;
% Inverse transform and mask result
t2 = ifft2((t1)); 
dZ = real(t2);
 dZ = dZ(sr+1:2*sr,sc+1:2*sc);   
dZ = (dZ).*dMask;
dZ(dMask==0)=NaN;

