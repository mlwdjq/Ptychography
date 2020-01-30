function [ud,fx,fy]=AngSpec(u,lambda,distance,pixsize)
%u为输入复振幅
%lambda为波长（单位：mm）
%distance为传播距离
%pixsize为采样间隔（单位：mm）
%ud为输出复振幅
%（fx,fy）为输出平面坐标
%%
[M,N]=size(u);
k=2*pi/lambda;
fxs=1/pixsize;
fys=1/pixsize;
[fx,fy]=meshgrid(linspace(-fys/2,fys/2-1/pixsize/N*(1-mod(N,2)),N),linspace(-fxs/2,fxs/2-1/pixsize/M*(1-mod(M,2)),M));
U=fftshift(fft2(u));
Ud=U.*exp(1i*k*distance*sqrt(1-(fx*lambda).^2-(fy*lambda).^2));
ud=ifft2(ifftshift(Ud));