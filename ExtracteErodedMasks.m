function erodedMasks = extractErodedMasks(cell_masks, cellIDs, erosionPixels)
    % Extract and erode individual cell masks
    % cell_masks: 2D label image from Cellpose (each cell has unique integer)
    % cellIDs: array of cell IDs to extract
    % erosionPixels: number of pixels to erode (1-2 recommended)
    
    erodedMasks = {};
    
    if isempty(cell_masks) || isempty(cellIDs)
        return;
    end
    
    se = strel('disk', erosionPixels);
    
    for i = 1:length(cellIDs)
        cellID = cellIDs(i);
        % Extract this cell's binary mask
        singleMask = (cell_masks == cellID);
        % Erode it to be slightly smaller
        erodedMask = imerode(singleMask, se);
        % Store it
        erodedMasks{i} = erodedMask;
    end
end
