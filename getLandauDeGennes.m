% assumes coordinate input of (NPART x NDIM) matrix
% second input is list of (contiguous) particles to be included in LDG 
function [S] = getLandauDeGennes(INPCRD,RANGE, PBC,LSIMBOX) %RANGE=1:nscaf for me
    input_coords = INPCRD(RANGE,:);
    sz = size(input_coords);
    npart = sz(1);
    ndim = sz(2);
    dir = zeros(npart-1, ndim);
    for i=1:npart-1
        dist = input_coords(i+1,:)-input_coords(i,:);
        if (PBC)
            dir(i,:) = dist - round(dist./LSIMBOX,0).*LSIMBOX;
        else 
            dir(i,:) = dist;
        end
    end
    QVal = zeros(npart-1,3,3);
    for k = 1:npart-1
        nsdkt = norm(squeeze(dir(k,:)));
        for i=1:3
            for j=1:3            
                dir1 = dir(k,i)/nsdkt;
                dir2 = dir(k,j)/nsdkt;
                QVal(k,i,j) = 3/2 * dir1*dir2 - 1/2 * double(i==j);
                %fprintf("Qval = %i\n",QVal);
            end
        end
    end
    Q = squeeze(mean(QVal));  
    eMax = max(eig(Q));
    S = eMax;
end