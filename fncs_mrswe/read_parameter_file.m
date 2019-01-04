function [vsParams] = read_parameter_file(parameterFile)
% -------------------------------------------------------------
% read_parameter_file
% Matlab Version:  R2017.b
%
% Purpose:
%   Read parameters from text file necessary for shear wave velocity
%   reconstruction.
%
% Inputs:
%   parameterFile - path to parameter file
%
% Outputs:
%   vsParams - returns parameters needed to calculate
%              shear wave velocity and the shear modulus.
%
%
% Lorne Hofstetter
% University of Utah
% ---------------------------------------------------------------

    % Parameters in text file
    paramNames = {'pathToData';     %Location of MR shear wave data
        'vsMax';                    %Upper shear wave velocity limit [m/s]
        'rho';                      %Density of sample [kg/m^3]
        'gaussFilterSize2D';        %Sice of gaussian filter [pixels]
        'gaussFilterSigma2D';       %Filter width [pixels]
        'numberWaveFronts';         %Number of shear wavefronts to detect
        'sliceIndices';             %Which slices for velocity calculation
    }';

    paramValues = {};               %Placeholder

 
    fid = fopen(parameterFile, 'r');

    count = 1;
    row  = fgetl(fid);
    while ischar(row)
        row = strip(row);       %Remove trailing and leading white space.
        
        %Skip comments
        if length(row) >= 2    
            if strcmp(row(1:2), '//')
                row = fgetl(fid); 
                continue
            end
        end
        
        if count <= 1
            paramValues{count} = row;
        else
            paramValues{count} = str2num(row);
        end
        
        count = count + 1;
        row = fgetl(fid);
        
    end
    
    vsParams = cell2struct(paramValues, paramNames ,2);
    fclose(fid);
end