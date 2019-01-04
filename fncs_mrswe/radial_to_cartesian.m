function [cartImg] = radial_to_cartesian(velItrans, img, xtheta, xradius, row, col)
% Function transforms image from radial to cartesian.
% origin of coordinate transformation is specified by xCenter and yCenter
% Inputs:
%	 -velItrans - this is the velocity map and xtheta and xradius specify position details
%    -img - this is ARFI image (only used for sizing purposes)
%    -xtheta - vector containing theta coords
%    -xradius - vector containing radius coords  
%    -row - speficies the central ARFI point row coord
%    -col - speficies the central ARFI point col coord
%
%
% Lorne Hofstetter
% October 20, 2017
% University of Utah

cartImg = zeros(size(img));

for r =1:size(cartImg,1)
    for c = 1:size(cartImg,2)
        h = c-col;
        w = r-row;
        tmpTheta = atan2(h,w)*(360/(2*pi));
        tmpRadius = sqrt(h^2 + w^2);
        [tmpValTheta, indxTheta] = min(abs(xtheta - tmpTheta));
        [tmpValRadius, indxRadius] = min(abs(xradius - tmpRadius));
        cartImg(r,c) = velItrans(indxTheta, indxRadius);
        
    end
end   