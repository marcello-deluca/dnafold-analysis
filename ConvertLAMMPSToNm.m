function [TRAJOUT] = ConvertLAMMPSToNm(TRAJIN, BOXSIZE)
    NSIM = length(TRAJIN);
    for k=1:NSIM
        SIM = TRAJIN{k};
        BOXSZ = BOXSIZE{k};
        
        if length(SIM)<length(BOXSZ)
            BOXSZ(end)=[];
        elseif length(BOXSZ)<length(SIM)
            SIM(end)=[];
        end
        LSIM = length(SIM);
        for i=1:LSIM
            SIM{i} = (SIM{i}-0.5).*BOXSZ(i);
        end
        TRAJIN{k}=SIM;
    end
    TRAJOUT = TRAJIN;
end