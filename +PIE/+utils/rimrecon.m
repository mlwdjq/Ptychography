
%   function out = rimrecon(xderiv, yderiv)
% xderiv assumes that xderiv = wav(x + s,y) - wav(x,y)
% and yderiv = wav(x,y+s) - wav(x,y).
% This means that the x shear is to the left, and the y shear is down. The
% rightmost values in the xderiv matrix are thus not used, as well as the
% topmost values in the yderiv matrix.

% i.e. the following matrix:
%
%      8     1     6
%      3     5     7
%      4     9     2
%
% would have xderiv =
%
%     -7     5     x
%      2     2     x
%      5    -7     x
%
% and yderiv =
%
%      x     x     x    revise to   x     x     x
%      5    -4    -1               -5     4     1
%     -1    -4     5                1     4    -5

% 
function [out, T, A] = rimrecon(xderiv, yderiv, unNormalizedAlphas, unNormalizedBetas, abc) %"IMPROVED BY GGG" VERSION
if nargin == 2
    unNormalizedAlphas = ones(size(xderiv));
    unNormalizedBetas = ones(size(yderiv));
    
end
if ~all(size(xderiv) == size(yderiv))
  error('X and Y derivatives must be the same size');
end

[sr,sc] = size(xderiv);


%% Modulo derivs 2pi

% xderiv = mod(xderiv + pi, 2*pi) - pi;%remove
% yderiv = mod(yderiv + pi, 2*pi) - pi;%remove

v=[xderiv(:);yderiv(:)];
%% Weight values
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



%% This code block populates the rimmer matrix T

T=improved_generateT(alphas, betas);
T(1,2:end) = 0;
T(1,1) = 1;

%% This code block populates the solution vector A

A=improved_generateA(alphas, betas);
% set first element to 0 for reference:
V = A*v;
V(1) = 0;%remove
%%

out = reshape(T\V,sr,sc);


function out=improved_generateT(al, be)
[n, m]=size(al);
[n2, m2] = size(be);

if (n~=n2) || (m~=m2)
    error('alpha and beta matrices must be the same size')
end

out=spalloc(n*m, n*m, 5*n*m);


al=-1*al(:);
out=out+diag(sparse(al(1:n*m-n)), -n);

be=sparse(be(:));
be(1:n:end) = 0;
out=out+diag(-be(2:end), -1);

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

al(n*m-n+1:end)=[];
f1=-al(:); %multiply in alpha
out= diag(sparse(f1), -n);

out=out-circshift(out, [n, 0]);
out = circshift(out, [-n,0]);
% f2=-al(:);
% c=zeros(n^2);
% c=c + diag(f2);
% c=circshift(c, [n, 0]);

%b
j1=sparse(-be(:));
u=spalloc(n*m, n*m, 2*n*m);
u=u+diag(j1);

betas=sparse(be(:));
b= diag(betas);
b=circshift(b,[1,0]);

b=b+diag(j1);

%b=[k b]

%b(:,n^2)=0;
b=circshift(b,[0,1]);
b(:,1:n:end)=0;
%k=zeros(n^2,1)

mt=[out b];

