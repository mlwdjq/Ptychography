% Sparse Rimmer reconstruction with weights removed temporarily (for
% simplicity)

function [out, T, A] = rimrecon_b(xderiv, yderiv) %"IMPROVED BY GGG" VERSION

if ~all(size(xderiv) == size(yderiv))
  error('X and Y derivatives must be the same size');
end

[sr,sc] = size(xderiv);

v=[xderiv(:);yderiv(:)];




%% This code block populates the rimmer matrix T

T = genT(sr, sc);
T(1,2:end) = 0;
T(1,1) = 1;

%% This code block populates the solution vector A

A = genA(sr, sc);
% set first element to 0 for reference:
V = A*v;
V(1) = 0;%remove
%%

out = reshape(T\V,sr,sc);


function T = genT(n, m)

N = n*m;


T=spalloc(n*m, n*m, 5*n*m);


al = -ones(N,1);
T=T+diag(sparse(al(1:n*m-n)), -n);


be = sparse(ones(N,1));
be(1:n:end) = 0;
T = T+diag(-be(2:end), -1);

T=T+T';

%diag
d=sum(T,2);
T=T-diag(d,0);


function mt=genA(n, m)

N = n*m;
al = ones(N - n,1);
f1=-al(:); %multiply in alpha
A2 = diag(sparse(f1), -n);
A2 = A2 - circshift(A2, [-n, 0]);

%b
dVec = ones(N,1);
dVec(n:n:end) = 0;
dVecm = ones(N,1);
dVecm(1:n:end) = 0;

A1 = diag(sparse(dVec)) - diag(sparse(dVecm(2:end)), -1);


mt= -[A2 A1];

