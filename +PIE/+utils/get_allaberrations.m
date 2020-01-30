%%波面生成
%版权所有 朱文华
function Iout=get_allaberrations(I,Piston,Tiltx,Tilty,Defocus,Astigx,Astigy,Comax,Comay,Sph)
[M,N]=size(I);
x=-1:2/(M-1):1;
y=-1:2/(N-1):1;
[X,Y]=meshgrid(y,x);
Iout=I+Piston+Comax*X.*(X.^2+Y.^2)+Comay*Y.*(X.^2+Y.^2)+Sph*(X.^2+Y.^2).^2+Astigx*(X.^2)+Astigy*(Y.^2)+Defocus*(X.^2+Y.^2)+Tiltx*X+Tilty*Y;

