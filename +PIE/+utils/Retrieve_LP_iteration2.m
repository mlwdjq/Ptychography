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

function dZ = Retrieve_LP_iteration2(dPhiX, dPhiY, dSx, dSy,dMask,iteration)

% Scale sheared wavefronts by shear maps to create proper difference
% functions
[sr,sc]=size(dPhiX);
s=1;
spx=dSx*(sc-1)/s;
spy=dSy*(sr-1)/s;
p = dPhiX/s;
q = dPhiY/s;

% Smooth boundary conditions by tiling flipped copies of sheared wavefronts
[hei, wid] = size(p);

p1 =  flipud(p);
p2 = -fliplr(p);
q1 = -flipud(q);
q2 = fliplr(q);

p = [-rot90(p,2) p1;p2 p];
q = [-rot90(q,2) q1;q2 q];

% Transform sheared wavefronts and divide by frequency variable

[h2,w2] = size(p);

dMasks=[dMask,dMask;dMask,dMask];

[U,V] = meshgrid(linspace(-1/2,1/2-1/w2,w2),linspace(-1/2,1/2-1/h2,h2));
dx=-1i*sin(2*pi*spx.*U);
dy=-1i*sin(2*pi*spy.*V);
LP0=abs(spx.*U)>=0.5|abs(spx.*V)>=0.5;
LP4=(lsianalyze.utils.topHatCosWin(2*sr,sr*0.5,sr*0.8));
LP3=abs(spx.*U)>=0.1|abs(spx.*V)>=0.1|(abs(dx)<1e-6&abs(dy)<1e-6);

% LP=mod(abs(dSx*(sc-1)/LP_width.*U),0.5)>=0.05&mod(abs(dSx*(sc-1)/LP_width.*U),0.5)<=0.45|mod(abs(dSx*(sc-1)/LP_width.*V),0.5)>=0.05&mod(abs(dSx*(sc-1)/LP_width.*V),0.5)<=0.45;
% LP=~LP;
% LP(abs(dSx*(sc-1)/LP_width.*U)<0.2&abs(dSx*(sc-1)/LP_width.*V)<0.2)=0;
% LP=LP|LP0;
% LP(h2/2+1,w2/2+1)=1;
% do iteration to remove edge effect
for j=1:iteration
    P = fftshift(fft2(p));
    Q = fftshift(fft2(q));
    Fxy=P+1i*Q;
%     Fxy=Fxy.*LP4;
    Fxy(LP0)=0;
    Wxy = ifft2(ifftshift(Fxy));
    Wx=real(Wxy);
    Wy=imag(Wxy);
    p(dMasks==0)=Wx(dMasks==0);
    q(dMasks==0)=Wy(dMasks==0);
end


t1 = (P+1i*Q)./(dx+1i*dy);

% remove spectrums outside the band limit
t1(isinf(t1)) = 0;
t1(isnan(t1))=0;

%  LP=((dSx*(sc-1).*U).^2+(dSy*(sr-1).*V).^2)>=0.25|(abs(dx)<1e-6&abs(dy)<1e-6);
LPx=abs(dSx*(sc-1).*U)>=0.5|(abs(dx)<1e-6&abs(dy)<1e-6);
LPy=abs(dSx*(sr-1).*V)>=0.5|(abs(dy)<1e-6&abs(dy)<1e-6);

% t1(LP0)=0;
t1=t1.*LP4;
% load('s=1.mat');
% t1(LP3)=t1s(LP3);
% Inverse transform and mask result
t2 = ifft2(ifftshift(t1));
dZ = real(t2);

dZ = dZ(hei+1:2*hei,wid+1:2*wid);
% dZ = (dZ).*dMask;
% dZ(dMask==0)=NaN;