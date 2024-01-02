function [TRAJOUT] = RemoveStaplesFromTrajectory(TRAJIN, NSCAF)
    NSIM = length(TRAJIN);
    for k=1:NSIM
        SIM = TRAJIN{k};
        sz = size(SIM);
        LSIM = sz(1);
        fprintf("%i/%i\n",k,NSIM)
        for i=1:LSIM
            FRAME = SIM{i};
            FRAME = FRAME(1:NSCAF,3:5);
            SIM{i} = FRAME;
            if mod(i,1)==0
                fprintf("%i\n",i);
            end
        end
        TRAJIN{k}=SIM;
    end
    TRAJOUT = TRAJIN;
end