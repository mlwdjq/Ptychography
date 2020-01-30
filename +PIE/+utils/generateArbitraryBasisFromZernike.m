%% this function is based on unit circle Zernike basis to generate an orthogonal basis that has arbitrary domain
% fitting operation: coef=pinv(basis)*wave; 
% transfer the coefficients to unit circle Zernike: coefZrn=(coef./delta)'*M
function [basis,basisInDomain,M,delta]=generateArbitraryBasisFromZernike(domain,x,y,dZrn)
% inital mask from domain
[sr,sc]=size(domain);
[th,r]=cart2pol(x,y);
mask=ones(sr,sc);
mask0=ones(sr,sc);
mask(domain==0|isnan(domain))=0;
mask0(r>1)=0;
delta=zeros(dZrn+1,1);
basisZ=zeros(sum(mask(:)),dZrn);
basisInDomain=zeros(sr,sc,dZrn);
C=zeros(dZrn); 

% generate unit circle Zernike basis
for i=1:dZrn+1
    afn = zgen([],i-1, 'fnr');
    temp=afn(r,th);
    delta(i)=sqrt(dot(temp(mask0==1),temp(mask0==1))/sum(mask0(:)));
    basisZ(:,i)=temp(mask==1)/delta(i);
end

% calculate the symmetric matrix
for i=1:dZrn+1
    for j=i:dZrn+1
        C(i,j)=dot(basisZ(:,i),basisZ(:,j))/sum(mask(:));
        if i<j
            C(j,i)=C(i,j);
        end
    end
end
% calculate conversion matrix
A=chol(C);
M=inv(A'); 

% generate the orthogonal basis
basis=(M*basisZ')';
temp=zeros(sr,sc);
for i=1:dZrn+1  
    temp(mask==1)=basis(:,i);
    basisInDomain(:,:,i)=temp;
end

