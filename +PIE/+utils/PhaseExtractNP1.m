% Using N+1 phase shifting interferograms to extract phase.
% reference: Design and assessment of symmetrical phase-shifting algorithms
function phase=PhaseExtractNP1(I)
N=length(I)-1;
alpha=zeros(size(I));
beta=zeros(size(I));
for n=0:N
    alpha(n+1)=sin(2*pi*n/N);
    if n>0&&n<N
        beta(n+1)=-cos(2*pi*n/N);
    else
        beta(n+1)=-1/2;
    end
end
if mod(N,2)==0
    k=1:N/2-1;
    c=1/N*sum((N-2*k).*sin(4*pi*k/N));
else
    k=1:(N-1)/2;
    c=1/N*sum((N-2*k).*sin(4*pi*k/N));
end
p=sum(alpha.*I+c*I(end)-c*I(1));
q=sum(beta.*I);
phase=atan2(p,q);
