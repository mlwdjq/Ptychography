%% this function is used to calculate overlap ratio of two shifted cicles 
% r is the radius of the circle
% d is the offset between two circles
function ratio = overlapRatio(r,d)
if d>=2*r
    ratio = 0;
    return;
end
s = sqrt(r.^2 - (d/2).^2); % half Chord length
thita = acos(d/2./r);
S1 = thita.*r.^2; % Sector area
S2 = d/2.*s; % triangle
ratio = (S1-S2)*2/pi./r.^2;