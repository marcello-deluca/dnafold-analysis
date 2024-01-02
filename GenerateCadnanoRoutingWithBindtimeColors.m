%function [] = GenerateCadnanoRoutingWithBindtimeColors(InputTimes,fileName)
    % NOTE: BE CAREFUL WITH INDICES
    InputTimes = mbt;
    
    %GenerateBindingTimeMap
    
    % Import and parse JSON file
    fileName = '4HB_end_staple_long.json'; % filename in JSON extension
    fid = fopen(fileName); % Opening the file
    raw = fread(fid,inf); % Reading the contents
    str = char(raw'); % Transformation
    fclose(fid); % Closing the file
    data = jsondecode(str); % Using the jsondecode function to parse JSON from string
    
    % Now process the data to create a plot showing the shape of the caDNAno
    % design
    
    % Find the 3' end first and then follow by 5' neighbors until end
    nrows = length(data.vstrands);
    for i = 1:nrows
        strand = data.vstrands(i);
        scaffolds = strand.scaf;
        staples = strand.stap;
        yLoc(i) = strand.row;
        HelixNumber(i) = strand.num;
        for j = 1:length(scaffolds)
             Neighbors = scaffolds(j,:);
    %         fprintf("For %i %i neighbors are %i %i %i %i\n", HelixNumber(i),j,Neighbors(1),Neighbors(2),Neighbors(3),Neighbors(4));
             FivePrimeNeigh = [Neighbors(1)+1,Neighbors(2)+1];
             ThreePrimeNeigh = [Neighbors(3)+1,Neighbors(4)+1];
    %         %fprintf("row %i. 3' neigh: %i. 5' neigh: %i.\n",i, ThreePrimeNeigh, FivePrimeNeigh);
            
            if isequal(ThreePrimeNeigh, [0 0]) && ~isequal(FivePrimeNeigh, [0 0])
               fprintf("found 3' end at %i %i\n",HelixNumber(i)+1,j);
               ThreePrimeEnd = [i j];
            end
        end
    end
    % Now Have Three Prime End, record it.
    ScaffoldRoute=[];
    ScaffoldRoute = [ScaffoldRoute; ThreePrimeEnd];
    
    % Current location is 3' end in internal notation
    CurrLoc = ThreePrimeEnd;
    % CurrentNeighbors = data.vstrands(CurrLoc(1)).scaf(j,:);
    CurrentNeighbors = data.vstrands(CurrLoc(1)).scaf(CurrLoc(2),:);
    FivePrimeNeighbor = CurrentNeighbors(1:2);
    FivePrimeNeighborInternal = computeFivePrimeNeighborInternal(FivePrimeNeighbor,data);
    %FivePrimeNeighborInternalIndices = computeFivePrimeNeighborInternal(FivePrimeNeighbor,data);
    
    %Follow the scaffold all the way to the 5' end
    iters = 0;
    while ~isequal(FivePrimeNeighbor,[-1,-1])
        FivePrimeNeighborInternal = computeFivePrimeNeighborInternal(FivePrimeNeighbor,data);
        CurrLoc = FivePrimeNeighborInternal;
        ScaffoldRoute = [ScaffoldRoute;FivePrimeNeighborInternal];
        %CurrLoc = FivePrimeNeighbor
        %Curr = data.vstrands(FivePrimeNeighbor(1)).scaf;
        %ScaffoldRoute = [ScaffoldRoute; FivePrimeNeighbor];
        CurrentNeighbors = data.vstrands(CurrLoc(1)).scaf(CurrLoc(2),:);
        FivePrimeNeighbor = CurrentNeighbors(1:2);
        
        
    %     iters = iters + 1;
    %     if mod(iters,100)==0
    %         fprintf ("iteration number %i\n",iters);
    %     end
    end
    
    x = [];
    y = [];
    z = [];
    
    % plotting part
    xLocations = ScaffoldRoute(:,2)  * 0.33; % 0.33nm/nt
    yLocations = max(ScaffoldRoute(:,1))-ScaffoldRoute(:,1);
    
    for i=1:length(ScaffoldRoute)
        fprintf("%i\n", data.vstrands.num);
        z(i) = 0;%data.vstrands(ScaffoldRoute(i,1)).row * 2.5; % 2.5nm / row
        y(i) = data.vstrands(ScaffoldRoute(i,1)).num * 2.5; %data.vstrands(ScaffoldRoute(i,1)).col * 2.5; % 2.5nm / column
    end
    y = transpose(y);
    zLocations = transpose(z);
    
    for i=1:length(InputTimes)
       MeanTimesPerNucleotide(8*(i-1)+1:8*(i-1)+8) = InputTimes(i) * .01E-9 * 10000 * 3314; 
    end
    x = xLocations;
    %y = yLocations;
    z = zLocations;
    col = transpose(MeanTimesPerNucleotide);
    patch([x;nan],[y;nan],[z;nan],[col;nan], 'edgecolor', 'interp'); 
    colorbar;
    colormap(jet);
    plt=Plot();
    
    axis equal
    xlabel("Horizontal Axis, cadnano");
    grid off
    axis off
    
    xlim([min(x)-10,max(x)+1]);
    plt.LineWidth=3;
%end

function [FivePrimeNeighborOut] = computeFivePrimeNeighborInternal(FivePrimeNeighborIn, data)
    assert(length(FivePrimeNeighborIn)==2);
    FivePrimeNeighborHelixInternalIndex = findStrandIndex(FivePrimeNeighborIn(1),data);
    InternalIndices=[FivePrimeNeighborHelixInternalIndex, FivePrimeNeighborIn(2)];
    scaffold = data.vstrands(InternalIndices(1)).scaf;
    Neighbors = scaffold(InternalIndices(2)+1,:);
    %FivePrime = Neighbors(1:2);
   % InternalIndices = findStrandIndex(Neighbors(1),data.vstrands.num);
    FivePrimeNeighborOut = InternalIndices+[0 1];
end

function [InternalStrandIndex] = findStrandIndex(InternalHelixNumberIn, data)
    l = length(data.vstrands);
    for i=1:l
       if data.vstrands(i).num==InternalHelixNumberIn
           fprintf ("strandnum %i, InternalHelixNumberIn %i\n",data.vstrands(i).num,InternalHelixNumberIn);
           InternalStrandIndex = i;
       end
    end
end
