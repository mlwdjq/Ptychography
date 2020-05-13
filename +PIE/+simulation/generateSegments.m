%% 10*10 segments
N = 40;
Ns = 10;
segs = cell(Ns);
for i= 1:Ns
    for j= 1:Ns
        I0 = zeros(N);
        I0((i-1)*4+1:i*4,(j-1)*4+1:j*4)=1;
        segs{i,j} = I0;
%         figure(2),imagesc(I0),drawnow;
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
save('../+utils/segs100_2.mat','segs');