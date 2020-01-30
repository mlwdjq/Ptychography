function p = elipticalHole(Dx,Dy, r, c)
D=max(Dx,Dy);
[X,Y] = meshgrid(  linspace(-1,1,D)  );
frame = zeros(D);
frame((D*X/Dx).^2+(D*Y/Dy).^2 <= 1) = 1;

p = pad2(frame, r,c);