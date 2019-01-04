function [peakIndices] = find_peaks(Itrans, pixelRange, nPeaks, signFirstPeak, debug)
% ========================================================
% Function finds the peakIndices. Required to calculate transient shear wave elastography
% Inputs:
	% - Itrans this is the transformed displacement image (or phase difference map) where
	%	row coordinates are angle in degrees and col coordinates are radial distance.
	% - pixelRange this is set by the maximum shear wave velocity (Number of pixels).
	% - nPeaks sets the number of displacement peaks [m/s]
	% - signFirstPeak speficies whether first peak is positive (1) or negative (0)
% Outputs:
	% - peakIndices: this is a "nPeaks" column vector specifying the detected peak indices
	%	these values are cricial to calculate shear wave velocity
% 
% Method: This algorithm works as follows. (1) It first extracts the optimal convolution 
% kernal from Itrans to identify the peaks. (2) Once the optimal convolution kernels has
% been obtained they are then used to the peak location for each radial ray. These peaks
% are then returned. The minVel, maxVel and nPeaks are requred for both parts (1) and (2).
% To improve robustness, noise in Itrans can be masked to prevent false positives, or 
% pixelRange may be sufficient to limit this noise

% As a final step polynomial fitting step to further refine localiztion of detected
% positions

% Lorne Hofstetter
% University of Utah
% =========================================================

    if nargin == 5
        debug = true;
    else
        debug = false;
    end
    

	%------------------- GENERATING SPECIFIC FILTERS FOR PEAK DETECTION ---------------------

	% Averaging all radial lines to get characteristic shapes
	ItransMean = mean(Itrans,1);
	peak = [];

	gaussFilter = fspecial('gaussian',[1 15],2); %setting filter of length 15 and sig=2 to detect shapes 
	
    if signFirstPeak == -1
		gaussFilter = -gaussFilter;
	end

	indxRange = 1:pixelRange;	%Setting search range for first peak
	
	% Finding all peaks in averaged curve
	for j=1:nPeaks
		peak(j) = locate_position(ItransMean, gaussFilter, indxRange);
		gaussFilter = -gaussFilter; %Alternating sign of gauss filter 
		indxRange = peak(j):(peak(j)+pixelRange); %Advancing search window for next peak search
	end

	% Finding Average sepation between peaks
	avePeakSeparation = floor(mean(diff(peak)));

	if(mod(avePeakSeparation,2))
		filterRange = avePeakSeparation +1;
	else
		filterRange = avePeakSeparation;
	end

	filterHalfRange = filterRange/2;

	% Generating the filters
	padSize = 20; %Pad before extracting filters to ensure sufficient range if extending before or after data
	padItransMean = padarray(ItransMean, [0 padSize]);

	filters = {};
	for j=1:nPeaks
		filters{j} = padItransMean(peak(j)+padSize-filterHalfRange:peak(j)+padSize+filterHalfRange);
    end

    if debug == true
        % Plotting Filters
        figure;
        for j=1:nPeaks
            subplot(nPeaks,1,j); plot(filters{j}); ylabel(num2str(j)); title('Filters');
        end
    end

	%------------------------- End of Filter Generation -----------------------

	% ------------------ Now Finding Peaks using correlation for Each Radial Projection ---------------------

	peakIndices = zeros(size(Itrans,1), nPeaks); % Empty array to store indices
	indxRange = 1:pixelRange;	%Setting search range for first peak
	
	% Finding peaks
	for j=1:nPeaks
		peakIndices(:,j) = locate_position(Itrans, filters{j}, indxRange);
		indxStart = floor(mean(squeeze(peakIndices(:,j))));
		indxRange = indxStart:(indxStart+pixelRange);
	end

	% --------------------------------------------------------------------------

	% ----  Now refining index location using polynomial fitting to find sub pixel location
	for k = 1:size(peakIndices,1)
		for j = 1:size(peakIndices,2)
			indexs = (peakIndices(k,j)-filterHalfRange):(peakIndices(k,j)+filterHalfRange);
			minVal = indexs(1);
			maxVal = indexs(end);

            % Checks to ensure no negative indexes
            indexs(indexs <1) = 1;
            
			CurveValues = Itrans(k,indexs);
            
            [p, ~, mu] = polyfit(indexs,CurveValues,3);

			% generating finely sampled function
			indexsFine = linspace(minVal, maxVal, 100);
            curveFine = polyval(p,indexsFine, [], mu);
            
			if mod(j+signFirstPeak,2)
				[val indxFind] = max(curveFine);
			else
				[val indxFind] = min(curveFine);
			end
			peakIndices(k,j) = indexsFine(indxFind);

		end
    end
    
    % Adding in center ARF index location - position is 0
    peakIndices = cat(2, zeros(size(peakIndices,1),1), peakIndices);
end



function [peakIndx] = locate_position(Itrans, Filter, indxRange)
% This function used by find_peaks to perform normalized correlation for peak detection.
% Inputs:
%	- Itrans: the transformed displacement image. 
%   - Filter: is the desired filter.
%   - indxRange: Specifies the index range to do final search for maximum in. 
% Outputs:
%   - Returns the location of the desired peak

	ItransPeak = normxcorr2(Filter, Itrans);

	% Cropping image to original size after the xcorr2
	[rSize, cSize] = size(Filter);
	rowCrop = (rSize-1)/2;
	colCrop = (cSize-1)/2;

	ItransPeak = ItransPeak(rowCrop+1:end-rowCrop, colCrop+1:end-colCrop);

	%Visualizing ItransPeak
	%figure; subplot(2,1,1); plot(sum(ItransPeak,1)); title('Peak Mean');
	%subplot(2,1,2); plot(ItransPeak(1,:)); title('theta = 180');

	% Finding the location of the peak
	[tmpPeakVal, peakIndx] = max(ItransPeak(:,indxRange),[],2);
	peakIndx = peakIndx + indxRange(1)-1;

end

