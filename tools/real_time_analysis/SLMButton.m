assignin("base","target_cells",target_cells);
writematrix(target_cells,'stimulationTargets.csv');

% NEW: Extract eroded masks for stimulation (add at the beginning)
erosionPixels = 2;  % Adjust as needed (1-2 pixels recommended)
targetCellMasks = extractErodedMasks(app.cell_masks, cellIDs, erosionPixels);
targetCellMasks_group2 = extractErodedMasks(app.cell_masks, cellIDs_group2, erosionPixels);

if app.SLM_ON && app.SLM == "BNS" % Boulder Nonlinear Systems
    if stimParams.seqMode
        if size(target_cells_group2,1) > 0 % if there are two groups for stimulation

            nTargets = 2; % hard coded for now, assuming there are only two groups

            % MODIFIED: Pass masks and exptVars
            [phaseMask1,msk1,I1] = stimulateBNS_SLM(target_cells, stimParams, targetCellMasks, app.exptVars);
            pause(1)
            [phaseMask2,~,~] = stimulateBNS_SLM(target_cells_group2, stimParams, targetCellMasks_group2, app.exptVars);
            pause(1)

            image_size = size(phaseMask1);

            images_flat(:,1) = reshape(phaseMask1,image_size(1)*image_size(2),1);
            images_flat(:,2) = reshape(phaseMask2,image_size(1)*image_size(2),1);
            
            image_sequence = reshape(images_flat,image_size(1)*image_size(2)*2,1);

            disp('Downloading images... please wait.');
            max_OD_frames = 20;
            download_error = calllib(app.dll_name, 'Download_image_sequence', app.board_number, nTargets, image_sequence, max_OD_frames);
            if download_error
                disp('Failed to download images.');
            else
                disp('Done downloading images.');
                app.StimulateButton.Enable = 1;
                imshow(msk1,'Parent', app.UIAxes8);
                imshow(I1,'Parent', app.UIAxes8_2);
            end

        else
        
            % MODIFIED: Pass individual cell masks
            [phaseMask1,msk1,I1] = stimulateBNS_SLM(target_cells(1,:), stimParams, {targetCellMasks{1}}, app.exptVars);
            pause(1)
            [phaseMask2,~,~] = stimulateBNS_SLM(target_cells(2,:), stimParams, {targetCellMasks{2}}, app.exptVars);
            pause(1)
            if nTargets > 2
                [phaseMask3,~,~] = stimulateBNS_SLM(target_cells(3,:), stimParams, {targetCellMasks{3}}, app.exptVars);
                if nTargets == 4
                    [phaseMask4,~,~] = stimulateBNS_SLM(target_cells(4,:), stimParams, {targetCellMasks{4}}, app.exptVars);
                end
            end

            image_size = size(phaseMask1);

            images_flat(:,1) = reshape(phaseMask1,image_size(1)*image_size(2),1);
            images_flat(:,2) = reshape(phaseMask2,image_size(1)*image_size(2),1);

            if nTargets == 4
                images_flat(:,3) = reshape(phaseMask3,image_size(1)*image_size(2),1);
                images_flat(:,4) = reshape(phaseMask4,image_size(1)*image_size(2),1);
                image_sequence = reshape(images_flat,image_size(1)*image_size(2)*4,1);
            elseif nTargets == 3
                images_flat(:,3) = reshape(phaseMask3,image_size(1)*image_size(2),1);
                image_sequence = reshape(images_flat,image_size(1)*image_size(2)*3,1);
            else
                image_sequence = reshape(images_flat,image_size(1)*image_size(2)*2,1);
            end

            disp('Downloading images... please wait.');
            max_OD_frames = 20;
            download_error = calllib(app.dll_name, 'Download_image_sequence', app.board_number, nTargets, image_sequence, max_OD_frames);
            if download_error
                disp('Failed to download images.');
            else
                disp('Done downloading images.');
                app.StimulateButton.Enable = 1;
                imshow(msk1,'Parent', app.UIAxes8);
                imshow(I1,'Parent', app.UIAxes8_2);
            end
        end
    end
end
