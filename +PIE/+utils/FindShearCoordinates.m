function [ShearX,ShearY]=FindShearCoordinates(lambda,z1,z2,T,halfD)
x(1)=0;
i=1;
while 1
i=i+1;
th=atan(x(i-1)/(z1+z2));
ph=asin(sin(th)+lambda/T);
if z1*tan(th)+z2*tan(ph)<=halfD
x(i)=z1*tan(th)+z2*tan(ph);
else 
    break;
end
end
xn(1)=0;
i=1;
while 1
    i=i+1;
th=atan(xn(i-1)/(z1+z2));
ph=asin(sin(th)-lambda/T);
if z1*tan(th)+z2*tan(ph)>=-halfD
xn(i)=z1*tan(th)+z2*tan(ph);
else 
    break;
end
end
xs=[fliplr(xn) x(2:end)];
ShearX=repmat(xs,length(xs),1);
ShearY=ShearX';