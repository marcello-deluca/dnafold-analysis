function [times] = ObtainTimeResolvedBindingStatuses(FILENAME,NSIM,NPART)
    times = zeros(NSIM,NPART);
    for i = 1:NSIM
        dirname = sprintf("%i",i);
        cd (dirname);
        DAT = importdata(FILENAME,' ');
        DAT(end,:)=[];
        t = DAT(:,1);
        NBOUND = DAT(:,2);
        PARTICLENUM = DAT(:,3);
        times(i,:) = OneSimulationTimeResolvedBindingStatuses(t,NBOUND, PARTICLENUM, NPART);
        cd ..
    end
end
