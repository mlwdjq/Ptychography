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

function dZ = Retrieve_LP_iteration1(dPhiX, dPhiY, dSx, dSy,dMask,iteration)

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

dx=-1i*sin(2*pi*dSx*(sc-1).*U);%(exp(2*1i*pi*dSx*(sc-1).*U)-exp(-2*1i*pi*dSx*(sc-1).*U))/2;
dy=-1i*sin(2*pi*dSy*(sr-1).*V);%(exp(2*1i*pi*dSy*(sr-1).*V)-exp(-2*1i*pi*dSy*(sr-1).*V))/2;
t1 = (P+1i*Q)./(dx+1i*dy);

% remove spectrums outside the band limit
t1(isinf(t1)) = 0;
t1(isnan(t1))=0;
LP=(abs(dx)<1e-6&abs(dy)<1e-6);
  LP=abs(dSx*(sc-1).*U)>=0.5|abs(dSy*(sr-1).*V)>=0.5|(abs(dx)<1e-6&abs(dy)<1e-6);
    LP3=abs(dSx*(sc-1).*U)>=0.4|abs(dSy*(sr-1).*V)>=0.4;
% LP=((dSx*(sc-1).*U).^2+(dSy*(sr-1).*V).^2)>=0.25|(abs(dx)<1e-6&abs(dy)<1e-6);
LPx=abs(dSx*(sc-1).*U)>=0.5|(abs(dx)<1e-6&abs(dy)<1e-6);
LPy=abs(dSx*(sr-1).*V)>=0.5|(abs(dy)<1e-6&abs(dy)<1e-6);
LP5=abs(dx+1i*dy);
LP5=LP5/max(LP5(:));
LP5(LP5>0.1)=1;
% mesh(LP5)
% t1(LP)=0;

LP2=mod(abs(dSx*(sc-1).*U),0.5)>=0.05&mod(abs(dSx*(sc-1).*U),0.5)<=0.45|mod(abs(dSy*(sr-1).*V),0.5)>=0.05&mod(abs(dSy*(sr-1).*V),0.5)<=0.45;
LP2=~LP2;
LP2(abs(dSx*(sc-1).*U)<0.1&abs(dSy*(sr-1).*V)<0.1)=0;


LP4=(lsianalyze.utils.topHatCosWin(2*sr,1.9/(dSx*(sc-1))*sr,1.9/(dSx*(sc-1))*sr));
LP6=zeros(2*sr);
s=lsianalyze.utils.topHatCosWin(2*sr,5,20);
for m=-3:4
    for n=-3:4
        if m~=0||n~=0
        LP6=LP6+circshift(s,[32*m,32*n]);
        end
    end
end
LP6=(1-LP6);
% Inverse transform and mask result
t2 = ifft2(ifftshift(t1));
dZ = real(t2);

dxs=-1i*sin(2*pi.*U);%(exp(2*1i*pi*dSx*(sc-1).*U)-exp(-2*1i*pi*dSx*(sc-1).*U))/2;
dys=-1i*sin(2*pi.*V);%(exp(2*1i*pi*dSy*(sr-1).*V)-exp(-2*1i*pi*dSy*(sr-1).*V))/2;
% dxs=dx;
% dys=dy;
% dxs(LP)=0;
% dys(LP)=0;
% do iteration to remove edge effect
dMasks=[dMask,dMask;dMask,dMask];
dZ0=dZ;
for j=1:iteration*0
    dWx=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dxs)));
    dWy=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dys)));
    p(dMasks==0)=dWx(dMasks==0)*dSx*(sc-1);
    q(dMasks==0)=dWy(dMasks==0)*dSx*(sc-1);
    P = fftshift(fft2(p/dSx/(sc-1)));
    Q = fftshift(fft2(q/dSx/(sc-1)));
    t1 = (P+1i*Q)./(dxs+1i*dys);
    t1(isinf(t1)) = 0;
    t1(isnan(t1))=0;
     t1(LP)=0;
    % Inverse transform and mask result
    t2 = ifft2(ifftshift(t1));
    dZ = real(t2);
%     dZ=dZ/dSx/(sc-1);
end
t3=t1;

dZ=dZ0;
for j=1:iteration
    dWx=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dx)));
    dWy=real(ifft2(ifftshift(fftshift(fft2(dZ)).*dy)));
    p(dMasks==0)=dWx(dMasks==0);
    q(dMasks==0)=dWy(dMasks==0);
    P = fftshift(fft2(p));
    Q = fftshift(fft2(q));
    t1 = (P+1i*Q)./(dx+1i*dy);
    t1(isinf(t1)) = 0;
    t1(isnan(t1))=0;
%       t1(LP2)=0;
t1=t1.*LP4;
    % Inverse transform and mask result
    t2 = ifft2(ifftshift(t1));
    dZ = real(t2);
%     dZ=dZ/dSx/(sc-1);
end

t4=t1;

t4(LP3)=t3(LP3);
t2 = ifft2(ifftshift(t4));
% dZ = real(t2);
% figure(6),mesh(LP2)
% dZ=real(dZ);
dZ = dZ(hei+1:2*hei,wid+1:2*wid);
dZ = (dZ).*dMask;
dZ(dMask==0)=NaN;