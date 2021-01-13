function phase=DFTforPhaseExtraction(I,PS,frequency,Orders)
if nargin<4
    Orders=[-2,-1,0,1,2];
end
PS = PS(:);
Orders = Orders(:);
L=length(Orders);
S=(1-L)/2:(L-1)/2;
[~,Is]=sort(abs(S));
Index=S(Is);
num=length(Orders);
H=zeros(num);
[Sr,Sc,M]=size(I);
Is=zeros(Sr*Sc,M);
for i=1:M
    temp=I(:,:,i);
    Is(:,i)=temp(:);
end

for i=1:num
    for j=1:num
        H(i,j)=sum(exp(-1i*(Orders(i)+Index(j))*PS));
    end
end

phase=zeros(Sr,Sc);


F=Is*exp(-1i*PS*Orders');
C=H\F.';

% phase(phase==0)=atan2(imag(C(2,:)),real(C(2,:)));% 1st order
phase(phase==0)=atan2(imag(C(2*frequency,:)),real(C(2*frequency,:)));


