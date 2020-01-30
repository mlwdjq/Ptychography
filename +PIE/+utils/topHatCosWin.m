%% generate top-hat + cos window 

function win=topHatCosWin(N,topHat,cutoff)
if nargin<1
    N=100;
    topHat=25;
    cutoff=50;
end
% if cutoff<topHat+10
%     cutoff=topHat+10;
% end
[x,y]=meshgrid(([1:N]-round((N+1)/2))/N*2);
[~,r]=cart2pol(x,y);
k=2*(cutoff/N-topHat/N);
if k>0
    win=1/2+cos(2*pi*(r-topHat/N)/k)/2;
    win((r-topHat/N)/k>0.5)=0;
else
    win=zeros(N);
end
win(r<topHat/N)=1;

% win=zeros(N); % top-hat filter
% win(abs(x)<topHat/N&abs(y)<topHat/N)=1;
