% PROCESS:
% IMPORT COORDS
% CONVERT TO NM
% UNWRAP
% RG

%Takes 1xLSIM cell containing NPART x 3 matrices of coordinates in
%nanometers. Outputs RG and plots shaded error bar with SEM.

function [RG] = getTimeResolvedRG(COORDSOUT)
    sz = size(COORDSOUT);
    nsim = sz(1);
    LSIM = sz(2);
    dt = 10000 * 10 *.005E-9 * 3314;
    skips = 1;
    RG = zeros(nsim,length(1:skips:LSIM));
    for simnum = 1:nsim
        for i=1:skips:LSIM                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            COORDS = COORDSOUT{simnum,i};
            COORDS = COORDS(:,1:3);
            RGCURR = computeRG(COORDS);
            RG(simnum,floor(i/skips)+1) = RGCURR;
            if mod(i,1000)==0
                fprintf("%f/%f\n",i,29178);
            end
        end
    end
    % figure();
    % hold on
    % for p = 1:nsim
    %     plot(S(p,:));
    % end
    t = (1:skips:LSIM) .* dt;
    MeanRG = mean(RG);
    SEM = std(RG)/sqrt(nsim);
    %hold on
   % stepsize = 1;
    %figure();
    %plot(t,RG)
    %shadedErrorBar(t(1:stepsize:end),MeanRG(1:stepsize:end),SEM(1:stepsize:end), 'lineprops', '-red');
   % ylabel("RG");
    %xlabel("time (s)");
end