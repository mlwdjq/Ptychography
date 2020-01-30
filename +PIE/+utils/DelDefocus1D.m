function [ wf1 , coef] = DelDefocus1D( wf0 ,x)
% remove defocus
    wf0=wf0(:);
    L=length(wf0);
    if nargin ==1
        x=linspace(-1,1,L);    
    end
    x=x(:);
    basis=[ones(L,1), x, x.^2];
    
    coef=pinv(basis)*wf0;
    wf1 = wf0 - coef(3)*x.^2-coef(2)*x- coef(1);
end
