function [PSx,PSy]=FindPhaseShift(I)
I1=I(:,:,1,1);
[M,N,P,Q]=size(I);
f=fftshift(fft2(I1));
fx=abs(f);
fy=fx;
delta=5;
fx(:,1:round(N+1)/2+delta)=0;
fy(1:round(M+1)/2+delta,:)=0;
[px1,px2]=find(fx==max(fx(:)));
[py1,py2]=find(fy==max(fy(:)));
PSx=zeros(P,Q);
PSy=zeros(P,Q);
for i=1:P
    for j=1:Q
        f=fftshift(fft2(I(:,:,i,j)));
        ph=atan2(imag(f),real(f));
        PSx(i,j)=ph(px1,px2);
        PSy(i,j)=ph(py1,py2);
    end
end