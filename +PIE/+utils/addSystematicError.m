function wInt = addSystematicError(wInt,dMaxPhoton,dShot2Shot,dcFlare)
% add dcFlare
wInt = (1-dcFlare)*wInt + max(wInt(:))*dcFlare;
% photon noise
if (dMaxPhoton > 0)
    wInt = wInt./max(wInt(:)) * dMaxPhoton;
    wInt = wInt + sqrt(wInt).*randn(size(wInt));
end
% Shot to shot
wInt = wInt * (1 + dShot2Shot/100 * randn(1));

% % fix negative value
% wInt(wInt<0)=-wInt(wInt<0);
