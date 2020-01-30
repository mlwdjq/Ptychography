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

function dZ = Retrieve1Dflip(dPhiX, dSx,  dMask)

% Scale sheared wavefronts by shear maps to create proper difference
% functions
sr=length(dPhiX);
p = dPhiX;
 p2 = -fliplr(p);
p = [p2 p];


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
% Q = fftshift(fft2(q));
w2= length(P);

U = linspace(0,1-1/w2,w2);

 dx=2*1i*sin(2*pi*dSx.*U);
t1 = P./dx;
k=find(abs(dx)<1e-10);
for j=1:length(k)
    if k(j)==1|| k(j)==w2
        t1(k(j))=0;
    else
        t1(k(j))=(t1(k(j)+1)+t1(k(j)-1))/2;
    end
end
% t1 = (P+1i*Q)./(U+1i*V)./2/pi/1i;
% t1(isinf(t1)) = 0;
% t1(isnan(t1))=0;
% t1(abs(real(t1))>1e6)=0;
% t1(abs(imag(t1))>1e6)=0;
% Inverse transform and mask result
t2 = ifft((t1)); 
dZ = real(t2);
 dZ = dZ(sr+1:2*sr);   
dZ = (dZ).*dMask;
dZ(dMask==0)=NaN;

