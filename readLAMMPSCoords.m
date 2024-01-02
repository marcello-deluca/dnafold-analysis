function [BOXSZ,COORDS] = readLAMMPSCoords(trajfilename, npart,stepsperframe)
    fid = fopen(trajfilename,'r');
    %first header
    CURR = textscan(fid,'%s',27,'HeaderLines',0,'CollectOutput',1);
    BOXSZ = 2 * str2double(CURR{1}{16}); % get box size first frame
    CURR = textscan(fid,'%.0f%.0f%f%f%f%f%f%f%f%f%f%f%f%f',npart,'HeaderLines',0, 'CollectOutput',1);
    i=1;
    %SystemPrompt = sprintf("wc -l %s",trajfilename);
    %[status, cmdout]= system(SystemPrompt);
    FRAMES = [];%cell();
    while (~isempty(CURR{1})) 
        if mod(i,stepsperframe)==0
            fprintf("%i\n",i);
            FRAMES = [FRAMES;CURR];
        end
        CURR = textscan(fid,'%s',27,'HeaderLines',0,'CollectOutput',1);
        if (isempty(CURR{1}))
            BOXSZ(end)=[];
            break
        end
        if mod(i,stepsperframe)==0
            BOXSZ = [BOXSZ;2*str2double(CURR{1}{16})];
        end
        CURR = textscan(fid,'%.0f%.0f%f%f%f%f%f%f%f%f%f%f%f%f',npart,'HeaderLines',0, 'CollectOutput',1);
        i=i+1;
    end
    COORDS=FRAMES;
end