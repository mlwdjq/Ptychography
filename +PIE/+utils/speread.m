function f=speread(filename)

fid=fopen(filename,'r');
fseek(fid,42,'bof');
N= fread(fid,[1,1],'uint16');
fseek(fid,656,'bof');
M= fread(fid,[1,1],'uint16');
fseek(fid,1446,'bof');
num= fread(fid,[1,1],'int32');
f=cell(1,num);
for i=1:num
fseek(fid,4100+(i-1)*M*N,'bof');
f{i}=fread(fid,[M,N],'uint16');
end
fclose(fid);