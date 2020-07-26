N = 2000;
NA = 0.0875;
T_um = 170;
f_um = 3000;
lambda_um= 13.5e-1;
domian_um= 600;
[xp,yp] = meshgrid(linspace(-domian_um/2,domian_um/2,N));
pupil = zeros(N);
pupil(xp.^2+yp.^2<=(f_um*NA).^2) = 1;
theta=6;
s_um = lambda_um*f_um/T_um;
on_axis = 0;
QWLSI = exp(-1i*pi/lambda_um/f_um.*((xp-f_um*tand(theta)*on_axis+s_um).^2+yp.^2))+...
    exp(-1i*pi/lambda_um/f_um.*((xp-f_um*tand(theta)*on_axis-s_um).^2+yp.^2))+...
    exp(-1i*pi/lambda_um/f_um.*((xp-f_um*tand(theta)*on_axis).^2+(yp+s_um).^2))+...
    exp(-1i*pi/lambda_um/f_um.*((xp-f_um*tand(theta)*on_axis).^2+(yp-s_um).^2));
E = 4*exp(1i*2*pi/lambda_um*sind(theta).*xp)+QWLSI;
I=abs(E).^2;
th=32;
I(I<th)=0;
I(I>th)=1;
I = I.*pupil;
figure(2),imshow(I,[])

