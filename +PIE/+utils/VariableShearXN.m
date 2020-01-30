function Shear=VariableShearX(lambda,z1,z2,NA,T,M,N)
% M=500;
% N=500;
% lambda=0.0135;
% z1=7;
% z2=21000;
% NA=0.5;
% T=0.234;
l=tan(asin(NA))*(z1+z2);
x=linspace(-l,l,M);
y=linspace(-l,l,N);
[x,y]=meshgrid(y,x);
th=atan(x/(z1+z2));
%thy=atan(y/(z1+z2));
ph=asin(sin(th)-lambda/T);
xs=z1*tan(th)+z2*tan(ph);
%thx=atan(xs/(z1+z2));
Shear=(x-xs)/l;
%mesh(sin(thx),sin(thy),Shear)
