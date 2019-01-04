function [] = mrswe(pathToParameterFile)
% -------------------------------------------------------------
% mrwse - MR shear wave elastography
% Method described here in following manuscript:
% https://doi.org/10.1002/mrm.27647 

% Matlab Version:  R2017.b
% Running on macOS 10.13.6
%
% Purpose:
%   Reconstruct shear wave velocity maps and elastograms
%   from multi-point MR ARFI dataset.
%
% Optional Inputs:
%   pathToParameterFile - path to parameter file
%
%
% Lorne Hofstetter
% University of Utah, 2018
% ---------------------------------------------------------------

% =========== Setting Path to Dependencies =============
addpath(genpath('./fncs_mrswe/'));

% ====================== Loading Data ==================

    % Specify and check existence of parameter file
    if nargin < 1
        pathToParameterFile = pick_file();
    elseif ~check_filename_exist(pathToParameterFile)
        disp('Choose another parameter file.');
        pathToParameterFile = pick_file();
    end

    % Load Parameter file
    vsParams = read_parameter_file(pathToParameterFile);
    
    % Load MR Data
    load(vsParams.pathToData);
    
    disp('Calculating shear wave speed maps ...');

% ========= Initialize Placeholder Variables ==================
    [rr cc ss npts] = size(pDiff);
    
    vsMaps = zeros(rr, cc, length(vsParams.sliceIndices), npts, vsParams.numberWaveFronts);
    distanceBetweenPeaks = deal(vsMaps); 
    % [row, column, # slices to calc vs, # ARF points, # shear wave fronts]
    
    % store location of the ARF impulses
    row = zeros(length(vsParams.sliceIndices), npts);
    col = zeros(length(vsParams.sliceIndices), npts);

    indicesWavefrontLocations = {};

% ============================= Defining Constants ==============================
%   rGridIncrement:                      [-] increment factor for transformation to cylindrical coordinates
%                                           rGridIncrement * distPerPixel equals length of each voxel in 
%                                           the radial dimension
%
%   deltaT:                              [mm] time in ms between center of Motion Encoding Gradient Lobes
%   distPerPixel:                        [ms] pixel width
%   maxNumberPixelsBetweenWavefronts:    [-] defined by range of shear wave velocities to calculate over
% ---------------------------------------------------------------------------

    rGridIncrement = sqrt(2);
    deltaT = shd.tMeg + shd.spacingBetweenMEG;
    distPerPixel = shd.readoutFOV/(shd.BaseResolution * shd.ZFIvect(1)); 
    distPerPixelR = distPerPixel*rGridIncrement;
    maxNumberPixelsBetweenWavefronts = ceil((vsParams.vsMax * deltaT) / (rGridIncrement * distPerPixel));
    
    is_pixel_square(shd); % Check to ensure pixels are square inplane. 
    
    tArrival = [];
    tArrivalStruct = {};
    
% ========================= Calculate Shear Wave Velocity ====================
    % Iterating over number of slices and number of ARF points
    
    for ns=1:length(vsParams.sliceIndices) 

        for np=1:npts
            
            % Find center location of ARF point
            [tmpRow, tmpCol] = locate_center(pDiff(:,:,vsParams.sliceIndices(ns),np), ...
                                vsParams.gaussFilterSize2D, ...
                                vsParams.gaussFilterSigma2D );

            % Storing location of detected ARF impulse;
            row(ns,np) = tmpRow;
            col(ns,np) = tmpCol;


            % Transform pDiff and magnitude images to cylindrical coord where origin is central ARFI point
            [pDiffCylindrical magnitudeCylindrical, xtheta, xradius] = cartesian_to_radial( ...
                                pDiff(:,:,vsParams.sliceIndices(ns),np), ...
                                magnitude(:,:,vsParams.sliceIndices(ns)), ...
                                tmpRow, ...
                                tmpCol);
          
            % Identify peak trough positions of encoded shear waveferonts
            wavefrontPositionIndices = find_peaks(pDiffCylindrical, ...
                                maxNumberPixelsBetweenWavefronts, ...
                                vsParams.numberWaveFronts, ...
                                -1);
               
            % Generate images of the distance between each shear wave front in cylindrical coord
            distanceBetweenPeaksCylindrical = find_distances(wavefrontPositionIndices, ...
                                pDiffCylindrical, ...
                                distPerPixelR);
                           
            % Transforming distance maps back to cartesian
            for nw = 1:vsParams.numberWaveFronts
                distanceBetweenPeaks(:,:,ns,np,nw) = radial_to_cartesian(distanceBetweenPeaksCylindrical(:,:,nw), ...
                                pDiff(:,:,vsParams.sliceIndices(ns),np), ...
                                xtheta, ...
                                xradius, ...
                                tmpRow, ...
                                tmpCol);

            end
        end

        % Calculating the tArrival for distance between center ARFI point and next encoded wavefront
        [tmpTArrival, tmpTArrivalStruct] = calc_tarrival(squeeze(distanceBetweenPeaks(:,:,ns,:,:)), deltaT);
        tArrival(ns) = tmpTArrival;
        tArrivalStruct{ns} = tmpTArrivalStruct;
        
        % Generate velocity map for each point using tArrival and deltaT
        vs = cat(5, distanceBetweenPeaks(:,:,:,:,1)/tArrival(ns), ...
              distanceBetweenPeaks(:,:,:,:,2:end)/deltaT);
        vsFlattened = sum(vs,5);

        %Generate composite shear wave velocity map using median combination
        vsMedian = combine_velocity_map_median(vsFlattened, 2);

    end


% ========= Generate Result Structure to Return =========

    result = struct('magnitude', magnitude(:,:,vsParams.sliceIndices));
    result.pDiff = pDiff(:,:,vsParams.sliceIndices,:);
    result.vsMedian = vsMedian;
    result.vs = vs;
    result.vsFlattened = vsFlattened;
    result.tArrival = tArrival;
    result.tArrivalStruct = tArrivalStruct;
    result.row = row;
    result.col = col;


% =========== Generating one figure as sanity check ==========

% By default plotting slice associated with index closest to six
desiredSlice = 6;
[~, plotIndex] = min(abs(vsParams.sliceIndices-desiredSlice)); 
generate_useful_plots(result, plotIndex);
    
% =========== Saving output to file =============
    save_data(pathToParameterFile, result, vsParams, shd);

    % Display structure of saved data
    disp(result);
end



% ======== File Handling Supporting Functions ==============
function [fileExist] = check_filename_exist(pathToParameterFile)
    fileExist = (exist(pathToParameterFile, 'file') == 2);
    disp('File not found...');
end

function [filePath] = pick_file()
    [fileName pathName] = uigetfile('./param_files/*.txt','Select Parameter File to Read');
    filePath = [pathName fileName]; 
end

function [] = is_pixel_square(shd)
    readoutResolution = shd.readoutFOV/(shd.BaseResolution * shd.ZFIvect(1));
    phaseResolution = shd.phaseFOV/(shd.PhaseRes * shd.ZFIvect(2));
    if readoutResolution ~= phaseResolution
        warning(['Pixels are not square in inplane dimenions.', ...
                'current implementation only accurate for square pixels.']);
    end
end

function [] = save_data(pathToParameterFile, result, vsParams, shd)
    tmpString = split(pathToParameterFile,'/');
    filename = tmpString{end}(1:end-4);
    folderPath = './results/';
    if (exist('./results/') ~= 7)
        mkdir(folderPath);
    end
    fullpath = [folderPath filename '.mat'];
    save(fullpath, 'result', 'vsParams', 'shd')
end

function [] = generate_useful_plots(result, sliceIndex)
    
    velocityRange = [0 4];

    % ---- Displaying each ARFI displacement image and each velocity ring
    allCatPdiff = [];
    allCatVelMap = [];
    for ss=1:size(result.pDiff,3)
        catPdiff = cat(2, result.pDiff(:,:,ss,1), result.pDiff(:,:,ss,2), result.pDiff(:,:,ss,3), result.pDiff(:,:,ss,4)); 
        catVelMap = cat(2, result.vsFlattened(:,:,ss,1), result.vsFlattened(:,:,ss,2), result.vsFlattened(:,:,ss,3), result.vsFlattened(:,:,ss,4));

        for j=2:4
            n = 4*(j-1);
            catPdiff = cat(1,catPdiff, cat(2,result.pDiff(:,:,ss,n+1), result.pDiff(:,:,ss,n+2), result.pDiff(:,:,ss,n+3), result.pDiff(:,:,ss,n+4))); 
            catVelMap = cat(1,catVelMap, cat(2,result.vsFlattened(:,:,ss,n+1), result.vsFlattened(:,:,ss,n+2), result.vsFlattened(:,:,ss,n+3), result.vsFlattened(:,:,ss,n+4)));
        end
        allCatPdiff(:,:,ss) = catPdiff;
        allCatVelMap(:,:,ss) = catVelMap;
    end
    
    
    h{2} = figure('pos', [31 215 1320 582]);
    subplot(1,2,1)
        imagesc(allCatPdiff(:,:,sliceIndex), [-0.5 0.5]); colormap(gca, 'gray'); axis image; colorbar;
        title('Phase Maps');
    subplot(1,2,2)
        imagesc(allCatVelMap(:,:,sliceIndex), velocityRange); colormap(gca, 'jet'); axis image; colorbar;
        title('Shear Wave Speed Maps [m/s]');

        
    % Plotting the magnitude image and composite shear wave speed map
    h{3} = figure('pos', [31 215 1320 582]);
    subplot(1,2,1)
        imagesc(result.magnitude(:,:,sliceIndex)); axis image; colormap(gca, 'gray'); colorbar;
        title('Magnitude Image');
    subplot(1,2,2)
        imagesc(result.vsMedian(:,:,sliceIndex), velocityRange); axis image; colormap(gca, 'jet'); colorbar;
        title('Composite Shear Wave Speed Map [m/s]');
        
    
end


    

    

