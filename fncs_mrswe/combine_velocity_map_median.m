function velocityMap = combine_velocity_map_median(velocityMap3D, nAveragesThresh)
% ===============================================================
% Function take a simple median filter of overlapping values
% Inputs:
%	- velocityMap3D: this is the multi-point ARFI velocity map
%	- nAveragesThresh: Default to 0, shows all data. If = 2 it means that
%     only pixels with 3 or more values will be combined ot output map
% Outputs:
%   - Returns a median combined velocity map\
%
% Lorne Hofstetter
% November 1, 2017
% April 19, 2018 Modifided to handle multiple slices.
%
% ===============================================================

	if nargin < 2
		nAveragesThresh = 0
	end
	
	[r c s pts] = size(velocityMap3D);
	velocityMap = zeros(r,c,s);
	for j = 1:r
		for k = 1:c
            for slc = 1:s
			tmpVect = squeeze(velocityMap3D(j,k,slc,:));
			tmpVect(tmpVect == 0) = [];
			if(length(tmpVect) > nAveragesThresh)	%requiring three or moreoverlapping datapoints
				velocityMap(j,k,slc) = median(tmpVect);
			end
		end
	end
end
