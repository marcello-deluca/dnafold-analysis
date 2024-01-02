function [bound,FirstBindTime] = OneSimulationTimeResolvedBindingStatuses(t,n,part, npart)
    FirstBindTime = zeros(npart,1);
    spf = 1;
    dwell_threshold = 100;
    lsim = t(end);
    bound = false(npart,round(lsim/spf,0));
    for i=1:length(t)
        if i==1
          bound(part(i)+1,1:round(t(i)/spf,0))=false;
          bound(part(i)+1,round(t(i)/spf,0):end)=true;
        else
           if n(i)>n(i-1)
               %fprintf("%i\n",i);
               bound(part(i)+1,round(t(i)/spf,0):end)=true;
           elseif n(i)<n(i-1)
               bound(part(i)+1,round(t(i)/spf,0):end)=false;
           end
        end
        if mod(i,1000)==0
           fprintf("t %i/%i\n",i,length(t)); 
        end
    end
    %nplots = 10;
    %n_per_plot =  floor(npart / nplots);
    % for i=1:nplots
    %     subplot(nplots,1,i);
    %     plot(transpose(bound((i-1)*n_per_plot+1:i*n_per_plot,:)))
    %     ylim([0 1.1]);
    % end
    sz = size(bound);
    %  Now quantify when each particle dwells.
    for i=1:npart
        foundBindTime=false;
        tbind = 1;
        while (~foundBindTime)
          while bound(i,tbind)==false
              tbind=tbind+1;
              if tbind > sz(2)
                  foundBindTime=true;
                  break
              end
          end
          if tbind+dwell_threshold-1 < sz(2)
              if bound(i,tbind:tbind+dwell_threshold-1) == ones(1,dwell_threshold)
                foundBindTime = true;
              else
                  tbind = tbind+1;
              end
          end
          if ~foundBindTime && tbind+dwell_threshold-1>sz(2)
              foundBindTime=true;
              FirstBindTime(i)=sz(2);
          end
          if (tbind)>=sz(2)
            foundBindTime = true;
          end
        
        end
        FirstBindTime(i)=tbind;
        fprintf("%i/%i\n",i,npart)
    end
end
