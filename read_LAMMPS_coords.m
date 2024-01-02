% READS OXDNA SIMULATION AND OUTPUTS TRAJECTORY COORDS FOR VISUALIZATION
function [BOXSZ,COORDS] = read_LAMMPS_coords(TRAJFILENAME, NPART, readics)
    BOXSZ=[];
    fid = fopen(TRAJFILENAME,'r');
    NLINES = NPART;
    i=1;
    for k=1:6
        line = fgetl(fid);
    end
    boxsz = textscan(fid,'%f %f',1,'CollectOutput',true);
    BOXSZ(1)=boxsz{1}(2)-boxsz{1}(1);
    for k=1:3 %2 extra lines and ITEM: ATOMS id type xs ys zs etc.
        line = fgetl(fid);
    end
    %CURR = textscan(fid, '%s %f %f %f', 3,'HeaderLines',0, 'CollectOutput',true);
    if readics
        CURR = textscan(fid, '%.0f%.0f%f%f%f%f%f%f%f%f%f%f%f%f', NLINES, 'CollectOutput',true);
    else 
        CURR = textscan(fid, '%.0f%.0f%f%f%f', NLINES, 'CollectOutput', true);
        %K = cell2mat(CURR);
    end
    i=2;
    FRAME = [];
    while(~isempty(CURR{1}))
        if mod(i,1000)==0
            %fprintf("%i\n",i)
        end
        if mod(i,1)==0
            FRAME = [FRAME,CURR];
        end     
        throwout = textscan(fid,'',3,'CollectOutput',false);

        %for k=1:3
        %    line = fgetl(fid);
        %end

        if line==-1
            break;
        end
        boxsz = textscan(line,'%f %f');
        BOXSZ = [BOXSZ, boxsz{2}-boxsz{1}];

        for k=1:3
            line = fgetl(fid);
        end
        %CURR = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 9, 'HeaderLines', 0, 'CollectOutput', 1);
        if readics
            CURR = textscan(fid, '%.0f%.0f%f%f%f%f%f%f%f%f%f%f%f%f', NLINES,'CollectOutput',true);
        else
            CURR = textscan(fid, '%.0f%.0f%f%f%f', NLINES, 'CollectOutput',true);
           % K=cell2mat(CURR);
        end  
        %FRAME = [FRAME,CURR];
        %CURR = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', numlines, 'HeaderLines', 2, 'CollectOutput', 1);
        i=i+1;
    end
    COORDS = FRAME;
    while length(COORDS)>length(BOXSZ)
        COORDS(end)=[];
    end
end