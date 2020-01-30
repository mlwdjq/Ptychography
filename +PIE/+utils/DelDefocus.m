function [ out ] = DelDefocus( in )
% 消倾斜+消离焦
    
[hei,wid] = size(in);
[x,y] = meshgrid(linspace(-1,1,wid),linspace(-1,1,hei));
vp = ~isnan(in);
inX = x(vp(:));
inY = y(vp(:));
inZ = in(vp(:));
pol1 = inX.^2;
pol2 = inY.^2;
pol3 = inX.*inY;
pol4 = 2*pol1 + 2*pol2 - 1;
% A是一个对称矩阵，为减少计算量先计算上三角部分然后复制
A = [sum(pol4.^2) sum(inX.*pol4) sum(inY.*pol4) sum(pol4);
    0 sum(inX.^2) sum(pol3) sum(inX);
    0 0 sum(inY.^2) sum(inY);
    0 0 0 length(inX)];
A = A+triu(A,1)';
B = [sum(pol4.*inZ);sum(inX.*inZ);sum(inY.*inZ);sum(inZ)];
coef = A\B;

out = in - coef(1)*(2*x.^2+2*y.^2-1) - coef(2)*x - coef(3)*y - coef(4);
    
end