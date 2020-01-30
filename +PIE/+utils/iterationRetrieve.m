%% this function apply iteration process to the FFT reconstruction method, which has better performance if there is an aperture
function dZ=iterationRetrieve(dPhiX, dPhiY, dSx, dSy,dMask,iterations)
if nargin<6
    iterations=50;
end

[sr,sc]=size(dPhiX);
if ~all(dMask(:)~=0)
    for loop=1:iterations
        dZ = -lsianalyze.utils.RetrieveWithoutMask(dPhiX, dPhiY, dSx, dSy);
        dWx=(sc-1)*dSx/2*(circshift(dZ,[0 1])-circshift(dZ,[0 -1]));
        dWy=(sr-1)*dSy/2*(circshift(dZ,[1 0])-circshift(dZ,[-1 0]));
        dPhiX(dMask==0)=dWx(dMask==0);
        dPhiY(dMask==0)=dWy(dMask==0);
        res=[dWx(dMask==0)-dPhiX(dMask==0);dWy(dMask==0)-dPhiX(dMask==0)];
        rms=std(res);
        if rms<1e-5
            break;
        end
    end
    dZ(dMask==0)=NaN;
else
    dZ = -lsianalyze.utils.Retrieve(dPhiX, dPhiY, dSx, dSy, dMask);
end


