function [ wf1 , coef] = removeSphericalPropagationPhase1D( wf0,lambda_nm,thita,coef1)
% remove defocus
    wf0=wf0(:);
    L=length(wf0);
    thita=thita(:);
    basis=[2 * pi /lambda_nm .* cos(thita),ones(L,1)];
    coef=pinv(basis)*wf0;
    if nargin == 4
        coef(1)=coef1;
    end
    wf1 = wf0 - basis*coef;
end
