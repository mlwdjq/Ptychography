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
%% save object
object = object_amp.*exp(1i*object_phase);
% object = fftshift(object);
save('../../data/object/MET5Grating.mat','object');