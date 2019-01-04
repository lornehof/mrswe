function [rCoord, cCoord] = locate_center(img, filtSize, sigmaa, debug)
% Function uses a 2D gaussian filter to identify the center ARFI point in image
% Inputs:
%	- img is a 2D ARFI displacement image
%   - filtSize x filtSize is the size of the gaussian shape
%	- sigmaa - specifies the sigma of the gaussian distribution.
%  Both filtSize and sigmaa are optional and default values are below. 
% Lorne Hofstetter
% October 20, 2017

    if nargin == 4
        debug = true;
    else 
        debug = false;
    end
    
	if nargin < 3
		sigmaa = 2.5; 
	end
	if nargin < 2
		filtSize = 15; 
	end
	
	% Defining 2D gaussian filter
	kernel = fspecial('gaussian', filtSize, sigmaa);

	imgCorr = normxcorr2(kernel, img);
    if debug == true
        figure, surf(imgCorr), shading flat
    end
    
	% Identifying center position
	[C, I] = max(imgCorr(:));
	[rCoord, cCoord] = ind2sub(size(imgCorr), I);

	% Adjusting for padding in xcorr2
	rCoord = rCoord - (filtSize-1)/2;
	cCoord = cCoord - (filtSize-1)/2;
	
	end
