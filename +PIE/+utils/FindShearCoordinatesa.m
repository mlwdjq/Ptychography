function [ShearX,ShearY]=FindShearCoordinatesa(lambda,z1,z2,T,halfD)
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
xa(1)=0;
i=1;
while 1
i=i+1;
ph=atan(xa(i-1)/z2);
th=asin(sin(ph)-lambda/T);
if (z1+z2)*tan(th)>=-halfD
xa(i)=(z1+z2)*tan(th);
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
xna(1)=0;
i=1;
while 1
i=i+1;
ph=atan(xna(i-1)/z2);
th=asin(sin(ph)+lambda/T);
if (z1+z2)*tan(th)<=halfD
xna(i)=(z1+z2)*tan(th);
else 
    break;
end
end
xs=[fliplr(xa) x(2:end)];
xns=[fliplr(xn) xna(2:end)];
xs=(xs+xns)/2;
ShearX=repmat(xs,length(xs),1);
ShearY=ShearX';