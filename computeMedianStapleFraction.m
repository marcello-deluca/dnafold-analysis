function medfrac = computeMedianStapleFraction(FILENAME, npart, NSIM, complements, strandnum)
    for k=1:NSIM
        fprintf("computing for simulation %i\n",k);
        FOLDERNAME = sprintf("%i",k);
        cd (FOLDERNAME);
        DAT = importdata(FILENAME,' ');
        DAT(end,:)=[];
        t = DAT(:,1);
        n = DAT(:,2);
        part = DAT(:,3);
        spf = 100E3;
        lsim = t(end);
        bound = false(npart,round(lsim/spf,0));
        for i=1:length(t)
            if i==1
              bound(part(i)+1,1:round(t(i)/spf,0)+1)=false;
              bound(part(i)+1,round(t(i)/spf,0)+1:end)=true;
            else
               if n(i)>n(i-1)
                   bound(part(i)+1,round(t(i)/spf,0)+1:end)=true;
               elseif n(i)<n(i-1)
                   bound(part(i)+1,round(t(i)/spf,0)+1:end)=false;
               end
            end
            if mod(i,1000)==0
               fprintf("computing bind status matrix %i/%i\n",round(t(i)/spf,0),round(lsim/spf,0)); 
            end
        end
        for i=1:224
            bound(complements.PairedStapleParticle(i)+1,:) = bound(complements.ScaffoldParticle(i)+1,:);
        end
        nstrands = max(strandnum);
        sz = size(bound);
        npart = sz(1);
        lsim = sz(2);
        for i=1:lsim
            if mod(i,1000)==0
                fprintf("computing med frac stap bind %i/%i\n",i,lsim);
            end
            for j=2:nstrands
                %fprintf("%i/%i\n",j,nstrands);
                indices=find(j==strandnum);
                if sum(bound(indices,i))>0
                    %fprintf("%i bound at time %i\n",j,i);
                    boundstaples(j,i)=1;
                    fractbound(j,i) = sum(bound(indices,i)) / length(indices); 
                else
                    boundstaples(j,i)=0;
                    fractbound(j,i)=0;
                end
            end
        end
        
        sz = size(fractbound);
        
        for i=1:sz(2)
            medianfractbound(i)=mean(nonzeros(fractbound(:,i)));
        end
        if k==1
            currlen = length(medianfractbound);
        else
            if length(medianfractbound)<currlen
                currlen = length(medianfractbound);

            end
        end
        medfrac(k,1:currlen)=medianfractbound(1:currlen);
        cd ../
    end
end