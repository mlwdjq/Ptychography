% Sparse Rimmer reconstruction with weights removed temporarily (for
% simplicity)

function [out, T, A] = rimrecon_bm(xderiv, yderiv, unNormalizedAlphas, unNormalizedBetas) %"IMPROVED BY GGG" VERSION

if ~all(size(xderiv) == size(yderiv))
  error('X and Y derivatives must be the same size');
end

[sr,sc] = size(xderiv);

v=[xderiv(:);yderiv(:)];

a = .25;
b = 20;
c = .2;


f = genRimWeightFn(a,b,c);

% unNormalizedAlphas = f(unNormalizedAlphas/(max(unNormalizedAlphas(:))));
% unNormalizedBetas = f(unNormalizedBetas/(max(unNormalizedBetas(:))));


%% This code block populates the rimmer matrix T

T = genT(unNormalizedAlphas, unNormalizedBetas);
T(1,2:end) = 0;
T(1,1) = 1;

%% This code block populates the solution vector A

A = genA(unNormalizedAlphas, unNormalizedBetas);
% set first element to 0 for reference:
V = A*v;
V(1) = 0;%remove
%%

out = reshape(T\V,sr,sc);


function T = genT(al, be)
[n,m]=size(al);
N = n*m;


T=spalloc(n*m, n*m, 5*n*m);


al = -1*al(:);
T=T+diag(sparse(al(1:n*m-n)), -n);


be = sparse(be(:));
be(1:n:end) = 0;
T = T+diag(-be(2:end), -1);

T=T+T';

%diag
d=sum(T,2);
T=T-diag(d,0);


function mt=genA(al,be)
[n,m]=size(al);
N = n*m;
al(n*m-n+1:end)=[];
f1=-al(:); %multiply in alpha
A2 = diag(sparse(f1), -n);
A2 = A2 - circshift(A2, [-n, 0]);

%b
dVec = be(:);
dVec(n:n:end) = 0;
dVecm = be(:);
dVecm(1:n:end) = 0;

A1 = diag(sparse(dVec)) - diag(sparse(dVecm(2:end)), -1);


mt= -[A2 A1];

