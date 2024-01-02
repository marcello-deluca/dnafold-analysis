% Takes string with filename (file containing 2 columns ASCII data, first column
% time, second column number bound particles) and number of simulations, and
% outputs OUTDATA, concatenated data with number of bound species WRT
% time. Assumes data is organized from current directory as
% {1..10}/filename.
function OUTDATA = gatherKineticsData(FILENAME, NSIM)
    OUTDATA = [];
    for i=1:NSIM
        FOLDERNAME = sprintf("%i",i);
        cd (FOLDERNAME);
        DAT = importdata(FILENAME,' ');
        NBOUND = DAT(:,2).*(DAT(:,2)<=512);
        SZ = length(NBOUND);
        if i==1
            MIN = SZ;
            OUTDATA = NBOUND;
        else
            if SZ < MIN
                MIN = SZ;
                OUTDATA(MIN+1:end,:)=[];
            end
            OUTDATA = [OUTDATA, NBOUND(1:MIN)];
        end
        fprintf("loaded file %i\n",i);
        cd ..
    end
end