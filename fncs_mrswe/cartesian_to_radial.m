function [Itrans, ItransMag, xtheta, xradius] = cartesian_to_radial(img, magimg, row, col)
% Function transforms image to from cartesian to radial. 
% origin of coordinate transformation is specified by xCenter and yCenter
% Inputs:
%	 -img - this is ARFI image
%    -magimg - this is the magnitude image
%    -row - speficies the central ARFI point row coord
%    -col - speficies the central ARFI point col coord
%
% Outputs: 
%	- Itrans - this is the transformed phase image
%   - magimg - this is the transformed magnitude image
%   - xtheta - this is the theta coords
%	- xradius - this is the radius subdifisions
%
% Lorne Hofstetter
% October 20, 2017
%

count = 1;
theta = 0;
radius = 0;
val = 0;
imag = 0; % keeping track of image magnitude

% Representing each pixel as three datapoints: radius, theta, magnitude
for r =1:size(img,1)
    for c = 1:size(img,2)
        h = c-col;
        w = r-row;
        theta(count) = atan2(h,w)*(360/(2*pi));
        radius(count) = sqrt(h^2 + w^2);
        val(count) = img(r,c);
        imag(count) = magimg(r,c);
        count = count + 1;
    end
end

% Defining the grid size for radial image
thetaInc = 1;
radiusInc = sqrt(2); % The radial distance sqrt(2) 
nradiusInc = ceil(max(size(img))/radiusInc); %Determining span of radial coords
R = nradiusInc*radiusInc + radiusInc/2; 
xtheta = -180:thetaInc:180;
xradius = radiusInc/2:radiusInc:R;
Itrans = zeros(length(xtheta), length(xradius));
ItransMag = Itrans;


% Stepping through all radii positions and colleecting all corresponding pixel, their magnitude and angle

for j = 1:size(Itrans,2)
    tmpRadius = radius;
    tmpTheta = theta;
    tmpVal = val;
    tmpImag = imag;
    
    % Selecting values within desired range
    tmpTheta(tmpRadius > (radiusInc/2 + (j-1)*radiusInc)) = [];
    tmpVal(tmpRadius > (radiusInc/2 + (j-1)*radiusInc)) = [];
    tmpImag(tmpRadius > (radiusInc/2 + (j-1)*radiusInc)) = [];
    tmpRadius(tmpRadius > (radiusInc/2 + (j-1)*radiusInc)) = [];

    tmpTheta(tmpRadius <= (-radiusInc/2 + (j-1)*radiusInc)) = [];
    tmpVal(tmpRadius <= (-radiusInc/2 + (j-1)*radiusInc)) = [];
    tmpImag(tmpRadius <= (-radiusInc/2 + (j-1)*radiusInc)) = [];
    tmpRadius(tmpRadius <= (-radiusInc/2 + (j-1)*radiusInc)) = [];
    
    if isempty(tmpRadius)
        % if tmpRadius is empty continue to next iteration.
        continue;
    end
    
    % Now stepping through each radii and averaging all theta values. Or if none assigning nearest value in theta
    for k=1:size(Itrans,1)
        tmptmpTheta = tmpTheta;
        tmptmpVal = tmpVal;
        tmptmpImag = tmpImag;
        tmptmpRadius = tmpRadius;

        tmptmpRadius(tmptmpTheta > (thetaInc/2 + xtheta(k))) = [];
        tmptmpVal(tmptmpTheta > (thetaInc/2 + xtheta(k))) = [];
        tmptmpImag(tmptmpTheta > (thetaInc/2 + xtheta(k))) = [];
        tmptmpTheta(tmptmpTheta > (thetaInc/2 + xtheta(k))) = [];
        
        tmptmpRadius(tmptmpTheta <= (-thetaInc/2 + xtheta(k))) = [];
        tmptmpVal(tmptmpTheta <= (-thetaInc/2 + xtheta(k))) = [];
        tmptmpImag(tmptmpTheta <= (-thetaInc/2 + xtheta(k))) = [];
        tmptmpTheta(tmptmpTheta <= (-thetaInc/2 + xtheta(k))) = [];
        
        if isempty(tmptmpTheta)
           [fval findx] = min(abs(xtheta(k)-tmpTheta));
           Itrans(k,j) = tmpVal(findx);
           ItransMag(k,j) = tmpImag(findx);
        else
           Itrans(k,j) = mean(tmptmpVal); 
           ItransMag(k,j) = mean(tmptmpImag);
        end
        
    end
end