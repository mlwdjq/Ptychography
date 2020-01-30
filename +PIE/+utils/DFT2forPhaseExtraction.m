function [PhX,PhY]=DFT2forPhaseExtraction(I,PSx,PSy,Orders)
if nargin<4
    Orders=[-2,-1,0,1,2];
end
[Ordersx,Ordersy]=meshgrid(Orders);
Ordersx=Ordersx(:);
Ordersy=Ordersy(:);
L=length(Orders);
S=(1-L)/2:(L-1)/2;
[~,Is]=sort(abs(S));
Index=S(Is);
[Indexx,Indexy]=meshgrid(Index);
Indexx=Indexx(:);
Indexy=Indexy(:);
PSx=PSx(:);
PSy=PSy(:);

num=length(Ordersx);
H=zeros(num);
[Sr,Sc,M,N]=size(I);
Is=zeros(Sr*Sc,M*N);
k=1;
for j=1:N
    for i=1:M
        temp=I(:,:,i,j);
        Is(:,k)=temp(:);
        k=k+1;
    end
end

for i=1:num
    for j=1:num
        H(i,j)=exp(-1i*(Ordersx(i)+Indexx(j))*PSx).'*exp(-1i*(Ordersy(i)+Indexy(j))*PSy);
    end
end

PhX=zeros(Sr,Sc);
PhY=zeros(Sr,Sc);

F=Is*(exp(-1i*PSx*Ordersx').*exp(-1i*PSy*Ordersy'));
C=H\F.';

PhX(PhX==0)=atan2(imag(C(L+1,:)),real(C(L+1,:)));
PhY(PhY==0)=atan2(imag(C(2,:)),real(C(2,:)));


