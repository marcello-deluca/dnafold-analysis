function indices = getStrandIndices(strandNums)
    nstrands = max(strandNums);
    for k=1:nstrands
        indices = find(k==strandNums);
        minIndex(k)=min(indices);
        maxIndex(k)=max(indices);
    end
end
