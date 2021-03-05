%% elliptical segments
N = 96;
Ns = 16;
segs = cell(Ns,1);
[x,y] = meshgrid(linspace(-2,2,N));
[th,r] = cart2pol(x,y);
Ellipticity = tan(asin(0.55/8))/tan(asin(0.55/4));
rs = [0.125,0.5625,1];
ths = linspace(-pi,pi,9);
for i= 1:length(rs)-1
    for j = 1:length(ths)-1
    I0 = zeros(N);
    rE = sqrt(x.^2+(y/Ellipticity).^2);
    if i==1&&j~=1
        I0(rE>=rs(i)&rE<=rs(i+1)&th>ths(j)&th<=ths(j+1))=1;
    elseif i~=1&&j==1
        I0(rE>rs(i)&rE<=rs(i+1)&th>=ths(j)&th<=ths(j+1))=1;
    elseif i==1&&j==1
        I0(rE>=rs(i)&rE<=rs(i+1)&th>=ths(j)&th<=ths(j+1))=1;
    else
        I0(rE>rs(i)&rE<=rs(i+1)&th>ths(j)&th<=ths(j+1))=1;
    end
    segs{(i-1)*(length(ths)-1)+j} = I0;
    figure(2),imagesc(I0),drawnow;pause(0.1)
    end
end

segs = segs(:)';
%% 10*10 segments
N = 40;
Ns = 10;
segs = cell(Ns);
for i= 1:Ns
    for j= 1:Ns
        I0 = zeros(N);
        I0((i-1)*4+1:i*4,(j-1)*4+1:j*4)=1;
        segs{i,j} = I0;
         figure(2),imagesc(I0),drawnow;
    end
end
segs = segs(:)';

%% 10*10 segments 2
N = 100;
Ns = 10;
segs = cell(Ns);
for i= 1:Ns
    for j= 1:Ns
        I0 = zeros(N);
        I0((i-1)*10+1:i*10,(j-1)*10+1:j*10)=1;
        segs{i,j} = I0;
         figure(2),imagesc(I0),drawnow;
    end
end
segs = segs(:)';

%% save data
save('../+utils/segs16e_2.mat','segs');