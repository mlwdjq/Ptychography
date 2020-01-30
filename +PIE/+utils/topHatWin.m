%% generate top-hat window 

function win=topHatWin(N,topHat)
if nargin<1
    N=100;
    topHat=25;
end
[x,y]=meshgrid(([1:N]-round((N+1)/2))/N*2);
[~,r]=cart2pol(x,y);
win=zeros(N);
win(r<topHat/N)=1;
