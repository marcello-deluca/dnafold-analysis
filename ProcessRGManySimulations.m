function [RG,TRAJOUT] = ProcessRGManySimulations(TRAJIN,BOX,NSCAF)
    TEMP = RemoveStaplesFromTrajectory(TRAJIN,NSCAF); %OK
    TEMP = ConvertLAMMPSToNm(TEMP,BOX);% OK
    TEMP = UnwrapManyTrajectoriesdnaBD(TEMP,BOX);
    TRAJOUT = TEMP;
    RG = getTimeResolvedRG(TRAJOUT);
end