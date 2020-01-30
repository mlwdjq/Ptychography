function zfn=generateZernikeFunction(zernCouples,N,flag)
if nargin<3
    flag=1;
end

switch flag
    case 1
        % Build zernike lambda using coef-mat method:
        tic
        coefMatSum = [];
        for k = 1:size(zernCouples, 1)
            coefMat = zernCouples(k,2) * zgen([], zernCouples(k,1), 'coef-mat');
            
            sz1 = size(coefMatSum);
            sz2 = size(coefMat);
            if isempty(coefMatSum)
                coefMatSum = coefMat;
            elseif (sz1(1) ~= sz2(1) || sz1(2) ~= sz2(2))
                % need to resize smaller matrix to fit bigger one:
                tempMat = zeros(max([sz1(1), sz2(1)]), max([sz1(2), sz2(2)]));
                tempMat(1:sz1(1), 1:sz1(2)) = coefMatSum;
                tempMat(1:sz2(1), 1:sz2(2)) = tempMat(1:sz2(1), 1:sz2(2)) + coefMat;
                coefMatSum = tempMat;
            else
                coefMatSum = coefMatSum + coefMat;
            end
            
        end
        % Next build up handles:
        
        fhR = cell(1,size(coefMatSum, 1));
        fhCos = cell(1, ceil(size(coefMatSum, 2)/2 + 1/2));
        fhSin = cell(1, ceil(size(coefMatSum, 2)/2 + 1/2));
        
        for k = 1:length(fhR)
            fhR{k} = @(r) r.^(k - 1);
        end
        for k = 1:length(fhCos)
            fhCos{k} = @(th) cos(k*th);
        end
        for k = 1:length(fhSin)
            fhSin{k} = @(th) sin(k*th);
        end
        zfn = @(r,th) zeros(N);
        
        for k = 1:size(coefMatSum, 1)
            for m = 1:size(coefMatSum, 2)
                if coefMatSum(k,m) == 0
                    continue
                end
                if m == 1
                    zfn = @(r, th) zfn(r,th) + coefMatSum(k,m) * fhR{k}(r);
                elseif mod(m,2) == 0
                    zfn = @(r, th) zfn(r,th) + coefMatSum(k,m) * fhR{k}(r).*fhCos{m/2}(th);
                else
                    zfn = @(r, th) zfn(r,th) - coefMatSum(k,m) * fhR{k}(r).*fhSin{(m - 1)/2}(th);
                end
            end
        end
        
        fprintf('Building zernike lambda using coef-mat took %s\n', s2f(toc));
        
        
    case 2
        % Build zernike handle using function method
        tic
        zfn = @(r,th) zeros(N);
        for k = 1:size(zernCouples, 1)
            afn = zgen([], zernCouples(k,1), 'fnr');
            zfn = @(r, th) zfn(r,th) + zernCouples(k,2) * afn(r, th);
        end
        fprintf('Building zernike lambda took %s\n', s2f(toc));
end