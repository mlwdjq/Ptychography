%%��������
%��Ȩ���� ���Ļ�
function Iout=get_allaberration(X,Y,Tiltx,Tilty,Defocus,Astigx,Astigy,Comax,Comay,Sph)
Iout=Comax*X.*(X.^2+Y.^2)+Comay*Y.*(X.^2+Y.^2)+Sph*(X.^2+Y.^2).^2+Astigx*(X.^2)+Astigy*(Y.^2)+Defocus*(X.^2+Y.^2)+Tiltx*X+Tilty*Y;

