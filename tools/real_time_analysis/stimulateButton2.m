if app.SelectforstimulationCheckBox.Value && app.SelectforstimulationCheckBox_2.Value
    cellIDs_selected = find(app.groupId == 1);
    cellIDs_correlated = app.topRankedCellIDs;
    cellIDs = union(cellIDs_correlated,cellIDs_selected);
elseif app.SelectforstimulationCheckBox_2.Value
    cellIDs = find(app.groupId == 1);
elseif app.SelectforstimulationCheckBox.Value
    cellIDs = app.topRankedCellIDs;
else
    return;
end

cellIDs_group2 = find(app.groupId == 2);

target_cells = horzcat(app.xcRaw(cellIDs),app.ycRaw(cellIDs));
target_cells_group2 = horzcat(app.xcRaw(cellIDs_group2),app.ycRaw(cellIDs_group2));

% NEW: Extract eroded masks for stimulation
erosionPixels = 2;  % Adjust as needed (1-2 pixels recommended)
targetCellMasks = extractErodedMasks(app.cell_masks, cellIDs, erosionPixels);
targetCellMasks_group2 = extractErodedMasks(app.cell_masks, cellIDs_group2, erosionPixels);

if size(target_cells_group2,1) > 0
    nTargets = 2;
else
    nTargets = size(target_cells,1);
end

stimParams = struct();
stimParams.seqMode = app.EnableSequenceModeCheckBox.Value;
stimParams.interval = app.IntervalmsEditField.Value;
stimParams.niter = app.IterationsEditField.Value;

if ~(stimParams.seqMode) && (app.SLM == "BNS")
    
    stimParams.VMode = app.VortexModeEditField.Value;
    stimParams.VLens = app.VirtualLensEditField.Value;
    stimParams.GS = app.GSIterationsEditField.Value;
    stimParams.yShift = app.XShiftEditField.Value;
    stimParams.xShift = app.YShiftEditField.Value;
    stimParams.dll_name = app.dll_name;
    stimParams.board_number = app.board_number;
    
    assignin("base","target_cells",target_cells);
    writematrix(target_cells,'stimulationTargets.csv');
end

if app.SLM_ON && app.SLM == "BNS"
    if stimParams.seqMode
        disp('seqMode');

        timerObj = timer('Name','MyTimer',               ...
           'Period',stimParams.interval/1000,            ... 
           'StartDelay',0,                 ... 
           'TasksToExecute',stimParams.niter,           ... 
           'ExecutionMode','fixedSpacing', ...
           'TimerFcn',@(~,~)stimulate_RT(nTargets,app.dll_name,app.board_number)); 

        start(timerObj);

    else
        zernike = readmatrix('zernike_coeffs.csv');
        xy = readmatrix('phaseMaskCenter.csv');

        stimParams.xShift = xy(2);
        stimParams.yShift = xy(1);
        
        % MODIFIED: Pass masks and exptVars to stimulateBNS_SLM
        [~,msk,I] = stimulateBNS_SLM(target_cells, stimParams, zernike, targetCellMasks, app.exptVars);
        
        for iter = 1:stimParams.niter
            calllib('PiUsb', 'piSetShutterState', 1, app.stimShutter);
            pause(0.03);
            calllib('PiUsb', 'piSetShutterState', 0, app.stimShutter);
            pause(stimParams.interval/1000);
        end

        imshow(msk,'Parent', app.UIAxes8);
        imshow(I,'Parent', app.UIAxes8_2);
    end

elseif app.SLM == "DMD"
    
    stimParams.seqMode = app.EnableSequenceModeCheckBox_2.Value;
    stimParams.testMode = app.EnableTestModeCheckBox.Value;
    stimParams.radius = app.RadiuspixelsEditField.Value;
    stimParams.pulseDuration= app.PulseDurationmsEditField.Value;
    stimParams.interval = app.RepeatEverymsEditField.Value;
    stimParams.niter = app.IterationsEditField_2.Value;
    stimParams.background = app.BackgroundIntensityEditField.Value;
    stimParams.stimGroup = app.StimGroupEditField.Value;
    stimParams.G1Size = app.Group1SizeEditField.Value;

    stimulateNikonDMD(target_cells,target_cells_group2,stimParams);

elseif app.SLM == "PrairieLink"
 
    stimParams.DStim = app.DStimulationDropDown.Value;
    stimParams.zOffset = app.ZoffsetEditField.Value;
    stimParams.isSpiral = 'False';
    if app.PatternDropDown.Value == "Spiral"
        stimParams.isSpiral = 'True';
    end
    paramsUsedInStim();
    stimulatePrairieLink(target_cells, stimParams, app.exptVars.dimX, app.exptVars.dimY);
    
else
    disp('Stimulation method is not available');
end
