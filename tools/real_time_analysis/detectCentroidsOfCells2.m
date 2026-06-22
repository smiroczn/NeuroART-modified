if inputParams.CFIND == 1 % Manual
    % stdImg = stdImg_window(regImg,20);
    % normStdIMG = (stdImg - repmat(min(stdImg(:)),size(stdImg))) ./ repmat(range(stdImg(:)),size(stdImg));
    figure; imagesc(normMeanIMG); colormap('gray'); axis('square'); clim([0 0.8]);
    selectText = ['Click on Neuron Centers...' newline 'Press Enter after all the cells are selected' newline 'Press delete if you want to deselect the last selected cell'];
    if verbose, disp (selectText),end
    [xc, yc] = getpts; %  manually select centers of the neurons
    cell_centroids(:,1) = yc;
    cell_centroids(:,2) = xc;
    clear xc yc
elseif inputParams.CFIND == 4 % From File
    [file,path] = uigetfile; % select the .mat file which contains the variable, "ptsIdx": column 2 --> x, column 3 --> y coordinates
    if file == 0
        fprintf('No file was selected ... \n');
        return;
    else
        load([path file],'cell_centroids');
    end
elseif inputParams.CFIND == 3 % Cite-on
    ups = 2.0; % upscaling factor
    th = 0.5;  % threshold
    [~,cmdOutput] = system(['activate cite-on & python citeon.py' ...
        ' -x ' num2str(exptVars.dimX) ...
        ' -y ' num2str(exptVars.dimY) ...
        ' -n ' num2str(length(batchSettings.frameBlock)) ...
        ' -p ' inputParams.imageFile ...
        ' -u ' num2str(ups) ...
        ' -t ' num2str(th)]);
    T = readtable('cell_coordinates.csv');
    clear cell_centroids;
    cell_centroids(:,1) = T.Var2; % yc
    cell_centroids(:,2) = T.Var1; % xc 

elseif inputParams.CFIND == 5 % cellPose
    % Cross-platform Python path detection
    if ispc  % Windows
        % Try common Windows Python locations
        possiblePaths = {
            'python',  % if Python is in PATH
            'C:\Python313\python.exe',
            'C:\Python312\python.exe',
            'C:\Users\%USERNAME%\AppData\Local\Programs\Python\Python313\python.exe',
            fullfile(getenv('LOCALAPPDATA'), 'Programs', 'Python', 'Python313', 'python.exe')
        };
        
        pythonPath = '';
        for i = 1:length(possiblePaths)
            [status, ~] = system([possiblePaths{i} ' --version']);
            if status == 0
                pythonPath = possiblePaths{i};
                break;
            end
        end
        
        if isempty(pythonPath)
            error('Python not found. Please install Python and ensure it is in your PATH.');
        end
        
    elseif ismac  % macOS
        pythonPath = '/Library/Frameworks/Python.framework/Versions/3.13/bin/python3';
        
    else  % Linux
        pythonPath = 'python3';
    end
    
    % Build and execute command
    cmd = sprintf('%s cellPose.py -n %d -p "%s"', ...
        pythonPath, ...
        length(batchSettings.frameBlock), ...
        inputParams.imageFile);
    
    [status, cmdOutput] = system(cmd);
    
    if status ~= 0
        error('Cellpose Python script failed with output: %s', cmdOutput);
    end
    
    % --- Load centroids from CSV ---
    if isfile('centroids.csv')
        T = readtable('centroids.csv');
        clear cell_centroids;
        cell_centroids(:,1) = T.Y; % yc
        cell_centroids(:,2) = T.X; % xc
    else
        error('Cellpose did not generate centroids.csv');
    end
    
    % --- Load full label masks ---
    if isfile('masks.tif')
        cell_masks = imread('masks.tif');
        fprintf('Loaded %d cell masks from Cellpose\n', max(cell_masks(:)));
    else
        warning('Cellpose mask file "masks.tif" not found. cell_masks set to [].');
        cell_masks = [];
    end
% elseif inputConfig.python.use_python && inputParams.CFIND == 2
%     cell_centroids = locateCentroids(inputConfig.python.complex_dynamics_path, inputConfig.python.env);

else
    if ~verbose
        [caimanTextOutput,cell_centroids,~,~,~] = evalc('CaImAn_CellFinder(regImg,caimanParams)'); % suppressing printed statements
    else
        [cell_centroids,~,~,~] = CaImAn_CellFinder(regImg,caimanParams); % CaImAn cell finder
    end
end





