function [phaseMask,msk,I] = stimulateBNS_SLM(targetCells, stimParams, zernike, cellMasks, exptVars)

% targetCells: x, y coordinates of the cells to be stimulated
% (first column -> x, second column -> y)
% zernike: Zernike polynomial coefficients
% cellMasks: cell array of binary masks (optional, uses circles if empty)
% exptVars: structure with dimX and dimY fields

% Created by Dulara De Zoysa, UMD on 5/13/2022
% Modified to support cell masks

% Handle optional parameters for backward compatibility
if nargin < 4
    cellMasks = {};
end
if nargin < 5
    exptVars = struct('dimX', 512, 'dimY', 512); % default values
end

  if ~exist('targetCells','var') || isempty(targetCells)
      error('targetCells must be provided and nonempty.');
  end

nTargets = size(targetCells,1); % Total number of cells to be stimulated
slmdim = 1536;
mult_factor = 1.25; % 1.5 for high pixel mode
zoom_factor = 3/4; % 1/6; % high pixel mode

xoff = [];
yoff = [];
vm = [];
           
for itr = 1:nTargets
    xoff = [xoff; targetCells(itr,1)]; 
    yoff = [yoff; targetCells(itr,2)];
    vm = [vm, stimParams.VMode]; % vortex modes
end
            
opts=dhot_opts();
opts.vortex_modes= vm; 
opts.niters = stimParams.GS;
opts.lens = -1*stimParams.VLens;
opts.slm_dims = [1536 1536]*mult_factor; % high pixel mode
            
% load the calibration file and transform 2P image coordinates to DHOT coordinates which are sent to the SLM 
P = PMTToDHOT(xoff,yoff); 

tstart = tic;

% After the check for cellMasks
if ~isempty(cellMasks) && itr <= length(cellMasks) && ~isempty(cellMasks{itr})
    fprintf('Using actual cell mask for cell %d\n', itr);  % ADD THIS
    cellMask = cellMasks{itr};
    % ... rest of code
else
    fprintf('Using circular disk for cell %d (no mask available)\n', itr);  % ADD THIS
    im = otslm.simple.aperture(sz, sz(1)/300, 'shape', 'circle');
end


When you stimulate, you should see in the command window:

Using actual cell mask for cell 1
Using actual cell mask for cell 2
...

if stimParams.VMode < 0 % consider negative values as disk shaped patterns, VMode=0 --> points, VMode>0 --> vortices
    c= Config.BNS(); % Calculations are done using EasySLM, based on otslm package
    sz= c.res;
    I = ones(c.res(1),c.res(2)); % uniform incident beam
    
    slm= DHOT1(c);
    
    % Create pattern for each cell
    for itr = 1:size(P,1)
        % Check if we have actual cell masks
        if ~isempty(cellMasks) && itr <= length(cellMasks) && ~isempty(cellMasks{itr})
            % Use the actual cell mask
            cellMask = cellMasks{itr};
            
            % Transform mask to SLM space
            % Get cell center in image space
            [rows, cols] = find(cellMask);
            if ~isempty(rows)
                % Create mask pattern in SLM space
                % Note: This is a simplified transformation
                % You may need to adjust based on your coordinate system
                maskSize = max(max(rows) - min(rows), max(cols) - min(cols));
                scaleFactor = (sz(1)/300) / (maskSize / 2); % scale to match expected size
                
                % Resize mask to appropriate size for SLM
                resizedMask = imresize(cellMask, scaleFactor, 'nearest');
                
                % Create aperture from mask
                im = double(resizedMask);
                im = im / max(im(:)); % normalize
                
                % Pad or crop to standard size
                targetSize = round(sz(1)/300) * 2;
                if size(im,1) > targetSize
                    im = imresize(im, [targetSize, targetSize], 'nearest');
                else
                    padSize = floor((targetSize - size(im,1))/2);
                    im = padarray(im, [padSize padSize], 0, 'both');
                    if size(im,1) > targetSize
                        im = im(1:targetSize, 1:targetSize);
                    end
                end
                
                % Ensure it's in the right format
                im = im(1:min(end,sz(1)), 1:min(end,sz(2)));
                padded_im = zeros(sz);
                padded_im(1:size(im,1), 1:size(im,2)) = im;
                im = padded_im;
            else
                % Fallback to circle if mask is empty
                im = otslm.simple.aperture(sz, sz(1)/300, 'shape', 'circle');
            end
        else
            % No mask available, use circular disk (original behavior)
            im = otslm.simple.aperture(sz, sz(1)/300, 'shape', 'circle');
        end
        
        % Create CGH pattern
        slmCGH = CGH(c, 'incident', I);
        pattern_for_cell = slmCGH.cgh(im, 'alpha', 0.5, 'use_gpu', false, 'N', 10);
        
        % Add to SLM
        slm.add_tweezer(P(itr,:), pattern_for_cell, 0);
    end

    pattern = slm.dhot('alpha',1);
    m = pattern;

else
    [m, ~]=dhot(P, opts, 0);
end

opts.zernike = zernike;
[~, IMG]=dhot(P, opts, 1);

t = toc(tstart);
fprintf('Time taken for dhot calculations before updating the display is %.2f seconds \n',t);
            
xspm = stimParams.xShift;
yspm = stimParams.yShift;
            
m = add_virtualLens(m, opts);
            
ref1 = round(slmdim*(mult_factor - 1)/2);
ref2 = round(slmdim*(mult_factor + 1)/2);
msk = m(ref1 + 1 + xspm: ref2 + xspm, ref1 + 1 + yspm: ref2 + yspm);
            
t = toc(tstart);
fprintf('Time taken until after updating the phase mask display in the GUI is %.2f seconds \n',t);

phaseMask = im2uint8(msk);
            
t = toc(tstart); 
fprintf('Total time taken for calculations before sending the phase mask to SLM is %.2f seconds \n',t);
            
low = 500;
high = prctile(IMG(:),99.9999);
if high < low
    low = 0;
end

IM = mat2gray(IMG,double([low,high]));
                
crop1 = mult_factor*slmdim*(0.5 - zoom_factor/2);
crop2 = mult_factor*slmdim*(0.5 + zoom_factor/2);

I = IM(crop1 + 1 : crop2, crop1 + 1 : crop2);

t = toc(tstart);
fprintf('Time taken until after updating the beamlet display in the GUI is %.2f seconds \n',t);

if ~(stimParams.seqMode)
    calllib(stimParams.dll_name, 'Write_overdrive_image', stimParams.board_number, phaseMask);
end

end
