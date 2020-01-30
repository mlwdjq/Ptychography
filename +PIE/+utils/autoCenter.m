%% this function is use to characterize the centriod 
function dCenter = autoCenter(ceInt)
[sr, sc] = size(ceInt);
I = zeros(size(ceInt{1,1}));
for i = 1:sr
    for j = 1:sc
        I = I+ceInt{i,j};
    end 
end
% level=graythresh(I); % find gray threshhold 

bw=imbinarize(I,'global'); % to bw 'adptive'
bw = bwareafilt(bw,1); % choose the largest region
C=regionprops(bw,'Centroid'); % get centroid
dCenter=round(C.Centroid);
