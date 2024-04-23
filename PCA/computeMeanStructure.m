%import coords first
% function inputs 
function MeanStruct = computeMeanStructure(input_coords)
    referenceFrame = -1;
    % End inputs
    sz = size(input_coords);
    npart = sz(1);
    lsim = sz(2);
    ndim = sz(3);
    if (referenceFrame==-1) 
        referenceFrame = ceil(rand*lsim); %pick random frame as reference
        fprintf("Using frame %i for mean structure computation\n",referenceFrame);
    end
    REFERENCE_FRAME_COORDS = squeeze(input_coords(:,referenceFrame,1:3));
    Centroid = mean(REFERENCE_FRAME_COORDS);
    REFERENCE_FRAME_COORDS = REFERENCE_FRAME_COORDS - Centroid;
    TimedStruct = zeros(lsim,npart,3);
    ExitLoop = false;
    iter = 1;
    while(~ExitLoop)
        fprintf("iteration %i\n",iter);
        if iter>1
            prev_mean_struct = MeanStruct;
            REFERENCE_FRAME_COORDS = MeanStruct;
        end
        for t=1:lsim
            CURRENT_FRAME_COORDS = squeeze(input_coords(:,t,1:3));
            % Transform current coords to best fit with reference coords
            [d,Z,tr] = procrustes(REFERENCE_FRAME_COORDS, CURRENT_FRAME_COORDS, 'scaling', false, 'reflection', false); 
            %plot3(Z(1:48,1),Z(1:48,2), Z(1:48,3));
    
            TimedStruct(t,:,:) = Z;
    
            if (mod(t,10000)==0)
               fprintf("%i/%i\n",t,lsim); 
            end
        end
        MeanStruct = squeeze(mean(TimedStruct));
        if (iter>1)
            d = procrustes(MeanStruct, prev_mean_struct, 'scaling', false, 'reflection', false);
            fprintf("d = %i\n",d);
            if d < 1E-6
               fprintf("Converged to mean structure\n");
               ExitLoop = true;
            end
        end
        iter = iter + 1;
    end
end

