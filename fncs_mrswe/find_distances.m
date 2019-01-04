function [displacementMaps] = find_distances(positionIndices, pDiffTrans, distancePerPixelR)
    % Lorne Hofstetter
    % University of Utah
    
    % Defining empty placeholders
    displacementMaps = zeros(size(pDiffTrans,1), size(pDiffTrans,2), size(positionIndices,2)-1);
    
    for j = 1:size(positionIndices,2)-1
        displacements = (positionIndices(:,j+1) - positionIndices(:,j))*distancePerPixelR;  
        
        % Adding measured displacements to displacement map
        for jRad = 1:size(displacementMaps,2)
            for kAngle = 1:size(displacementMaps,1)
                if jRad >= positionIndices(kAngle,j) && jRad <= positionIndices(kAngle,j+1)
                    displacementMaps(kAngle, jRad, j) = displacements(kAngle);
                end
            end
        end
        
        
    end
end

