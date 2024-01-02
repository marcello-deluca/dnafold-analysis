function [TRAJOUT] = RemoveStaplesFromTrajectory_singlesimulation(TRAJIN, NSCAF)
    LSIM = length(TRAJIN);
    for i=1:LSIM
        FRAME = TRAJIN{i};
        FRAME = FRAME(1:NSCAF,3:5);
        TRAJOUT{i} = FRAME;
        if mod(i,1)==0
            fprintf("%i\n",i);
        end
    end
end