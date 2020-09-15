%% this script is used to simulate and save object

N= 1000;

%% EUV mask
object_amp = ones(N);
object_phase = zeros(N);
sizeRange = 12:4:40;
ampRange = 0:0.2:0.8;
phaRange = 0:pi/6:1.5*pi;
[posx,posy]= meshgrid(25:50:975);
pos= [posx(:),posy(:)];
m=1;
for i =1:length(sizeRange)
    for j =1:length(ampRange)
        for k =1:length(phaRange)
            object_amp(pos(m,1)+[-sizeRange(i)/2:sizeRange(i)/2],pos(m,2)+[-sizeRange(i)/2:sizeRange(i)/2])=ampRange(j);
            
            object_phase(pos(m,1)+[-sizeRange(i)/2:sizeRange(i)/2],pos(m,2)+[-sizeRange(i)/2:sizeRange(i)/2])=phaRange(k) ;
            m = m+1;
        end
    end
end

figure,mesh(object_amp);
figure,mesh(object_phase);

%% EUV amp mask
object_amp = zeros(N);
object_phase = zeros(N);
sizeRange = 12:4:40;
ampRange = 1;
phaRange = 0;
[posx,posy]= meshgrid(25:50:975);
pos= [posx(:),posy(:)];
m=1;
for i =1:length(sizeRange)
    for j =1:length(ampRange)
        for k =1:length(phaRange)
            object_amp(pos(m,1)+[-sizeRange(i)/2:sizeRange(i)/2],pos(m,2)+[-sizeRange(i)/2:sizeRange(i)/2])=1;
            
%             object_phase(pos(m,1)+[-sizeRange(i)/2:sizeRange(i)/2],pos(m,2)+[-sizeRange(i)/2:sizeRange(i)/2])=phaRange(k) ;
            m = m+1;
        end
    end
end
 
figure,mesh(object_amp);
% figure,mesh(object_phase);

%% vaccum

object_amp = ones(N);
object_phase = zeros(N);

%% amplitude grating
[x,y] = meshgrid(linspace(-1.11036,1.11036,N));

object_amp =4+cos(2*pi*5*x)+cos(2*pi*5*y);
object_amp(object_amp<4) =0;
object_amp(object_amp>=4) =1;

object_phase = zeros(N);
figure,mesh(object_amp)
%% phase grating
[x,y] = meshgrid(linspace(-1,1,N));

object_phase =cos(2*pi*5*x)+cos(2*pi*5*y);

object_amp = ones(N);
figure,mesh(object_phase)

%% Four quadrant phase
object_amp = ones(N);

object_phase = zeros(N);
object_phase(1:N/2,N/2:end)=pi/2;
object_phase(N/2:end,N/2:end)=pi;
object_phase(N/2:end,1:N/2)=pi/2*3;

figure,mesh(object_phase)

%% object for WDD
n = 128; 
I1 =single(imread('pears.png'));
I1 = crop2(I1,n,n);
I2 =single(imread('onion.png'));
I2 = crop2(I2,n,n);
[x,y] = meshgrid(linspace(-1,1,n));
[x2,y2] = meshgrid(linspace(-1,1,64));
I1 = interp2(x,y,I1,x2,y2);
I2 = interp2(x,y,I2,x2,y2);

mask = zeros(64);
mask(17:48,17:48) =1;
% I1 = I1.*mask;
% I2 = I2.*mask;
I2 = mat2gray(I2)*0.8+0.2;
I1 =(mat2gray(I1)-0.5)*pi/2;
I1 = I1.*mask;
I2 = I2.*mask;
I1 = [I1,I1,I1;I1,I1,I1;I1,I1,I1];
I2 = [I2,I2,I2;I2,I2,I2;I2,I2,I2];
object_phase = crop2(I1,128,128);
object_amp = crop2(I2,128,128);
figure, imagesc(object_phase)

%% contact
n = 5;
ns = 118;
object_amp =(pad2(ones(n),ns,ns)+1)/2;
object_phase = pad2(ones(n),ns,ns)*pi*0.8;
% object_amp=circshift(object_amp,[2,2]);
% object_phase=circshift(object_phase,[2,2]);
figure(2),imagesc(object_amp)

%% two lines
n = 3;
ns = 118;
object_amp =pad2(ones(8*n,n),ns,ns);
object_amp = circshift(object_amp,[0,-3]) + circshift(object_amp,[0,3]);
object_phase = zeros(ns);
% object_amp=circshift(object_amp,[2,2]);
% object_phase=circshift(object_phase,[2,2]);
figure(2),imagesc(object_amp)
%% three lines
n = 2;
ns = 145;
object_amp =pad2(ones(8*n,n),ns,ns);
object_amp = object_amp+circshift(object_amp,[0,-4]) + circshift(object_amp,[0,4]);
object_phase = zeros(ns);
% object_amp=circshift(object_amp,[2,2]);
% object_phase=circshift(object_phase,[2,2]);
figure(2),imagesc(object_amp)

%% pinhole
n = 10;
ns = 96;
object_amp =pad2(pinhole(n),ns,ns);
object_phase = pad2(pinhole(n),ns,ns)*pi/4;
% object_amp=circshift(object_amp,[n/2,n/2]);
% object_phase=circshift(object_phase,[n/2,n/2]);
figure(2),imagesc(object_amp)

%% triangles
n = 15;
ns = 170;
s= ones(n);
[x,y]=meshgrid(linspace(-1,1,n));
s(x+2*y>1)=0;
s(x-2*y>1)=0;
object_amp =pad2(ones(n),ns,ns);
object_phase = pad2(ones(n),ns,ns)*pi/4;
figure(2),imagesc(s),colorbar

%% save object
object = object_amp.*exp(1i*object_phase);
% object = fftshift(object);
save('../../data/object/threeLine_4pixPitch_145.mat','object');