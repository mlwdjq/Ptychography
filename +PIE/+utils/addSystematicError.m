function wInt = addSystematicError(wInt,dMaxPhoton,dShot2Shot,dcFlare,segs)
% add dcFlare
wInt = (1-dcFlare)*wInt + max(wInt(:))*dcFlare;
% photon noise
if (dMaxPhoton > 0)
    wInt = wInt * dMaxPhoton;
    if ~isempty(segs)
            for k=1:length(segs)
                wInt = wInt + sqrt(wInt.*segs{k})* randn(1);
            end
    else
        wInt = wInt + sqrt(wInt).*randn(size(wInt));
    end
    wInt = wInt / dMaxPhoton;
end
% Shot to shot
wInt = wInt * (1 + dShot2Shot/100 * randn(1));

% % fix negative value
% wInt(wInt<0)=-wInt(wInt<0);
