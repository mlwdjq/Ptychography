

% 
function [out, T] = rimrecon3b(xderiv, yderiv, unNormalizedAlphas, unNormalizedBetas, abc) %"IMPROVED BY GGG" VERSION
if nargin == 2
    unNormalizedAlphas = ones(size(xderiv));
    unNormalizedBetas = ones(size(yderiv));
    
end
if ~all(size(xderiv) == size(yderiv))
  error('X and Y derivatives must be the same size');
end

[sr,sc] = size(xderiv);



v=[xderiv(:);yderiv(:)];

if nargin == 4
a = .25;
b = 20;
c = .2;
else
    a = abc(1);
    b = abc(2);
    c = abc(3);
end

f = genRimWeightFn(a,b,c);
alphas = f(unNormalizedAlphas/(max(unNormalizedAlphas(:))));
betas = f(unNormalizedBetas/(max(unNormalizedBetas(:))));




T=improved_generateT(alphas, betas);
T(1,2:end) = 0;

A=improved_generateA(alphas, betas);
V=A*v;
% set first element to 0 for reference:
V(1) = 0;%remove

out = reshape(T\V,sr,sc);

out = reshape(pinv(full(T))*V,sr,sc);


function out=improved_generateT(al, be)
[n, m]=size(al);
[n2, m2] = size(be);

if (n~=n2) || (m~=m2)
    error('alpha and beta matrices must be the same size')
end

out=spalloc(n*m, n*m, 5*n*m);

al=-1*al(:);
out=out+diag(sparse(al(1:n*m-2*n)), -2*n);

be=sparse(be(:));
be(n-1:n:end) = 0;
be(n:n:end) = 0;
out=out+diag(-be(1:end-2), -2);

out=out+out';

%diag
d=sum(out,2);
out=out-diag(d,0);


function mt=improved_generateA(al, be)
[n,m] = size(al);
[n2, m2] = size(be);

if (n~=n2) || (m~=m2)
    error('alpha and beta matrices must be the same size')
end

%a
f1=-al(:);
f1(n*m-2*n+1:end)=0;
al(n*m-2*n+1:end)=[];
f2=al(:); %multiply in alpha
out=spalloc(n*m, n*m, 2*n*m);
out=out+ diag(sparse(f2), -2*n)+diag(sparse(f1), 0);

% f2=-al(:);
% c=zeros(n^2);
% c=c + diag(f2);
% c=circshift(c, [n, 0]);

%b
be(n2-1:n2,:)=0;
f3=be(1:n2*m2-2);
j1=sparse(-be(:));
u=spalloc(n*m, n*m, 2*n*m);
u=u+diag(j1)+diag(sparse(f3), -2);

% betas=sparse(be(:));
% b= diag(betas);
% b=circshift(b,[1,0]);
% 
% b=b+diag(j1);
% 
% %b=[k b]
% 
% %b(:,n^2)=0;
% b=circshift(b,[0,1]);
% b(:,1:n:end)=0;
%k=zeros(n^2,1)

mt=[out u];

