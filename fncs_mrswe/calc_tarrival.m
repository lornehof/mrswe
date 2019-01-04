function [tarrival tArrivalStruct] = calc_tarrival(distanceBetweenPeaks, deltaT)
% puropse of this algorithm is to calculate the arrival time of the 
% initial shear wavefront from the formation of the ARF pulse to the 
% first shear wavefront. This tarrival is necessary to calculate velocity
% between the central point and the next encoded wavefront. 
%
% Lorne Hofstetter
% University of Utah, 2018

    velocityMaps = distanceBetweenPeaks(:,:,:,2:end)/deltaT;
    nRings = size(velocityMaps,4);

    
    % Generating a tArrival scatter for each point and each ring.
    tArrivalStruct = {};
    for j=1:nRings
        for ncenter = 1:size(velocityMaps,3)
            tArrivalTemp = [];
            for nouter = 1:size(velocityMaps,3)
                tArrivalTempInnerLoop = distanceBetweenPeaks(:,:,ncenter,1) ./ ...
                        velocityMaps(:,:,nouter,j);
                    
                % Removing divide by zero values.
                tArrivalTempInnerLoop(isnan(tArrivalTempInnerLoop))=[];
                tArrivalTempInnerLoop(isinf(tArrivalTempInnerLoop))=[];
                tArrivalTempInnerLoop(tArrivalTempInnerLoop == 0) = [];
                
                % Appending values
                tArrivalTemp = [tArrivalTemp tArrivalTempInnerLoop];
            end
            
            % Generating tArrival for each point using each ring
            tArrivalStruct{ncenter,j} = tArrivalTemp;
        end
    end
    
    [tarrivalMean, tArrivalMedian, tArrivalSigma] = calcAllData(tArrivalStruct); 
    
    % Usng median arrival time from all data points
    tarrival = tArrivalMedian;
end

function [tArrivalMean, tArrivalMedian, tArrivalSigma] = calcAllData(tArrivalStruct)
    tArrivalData = [];
    for j=1:size(tArrivalStruct,1)
        for k=1:size(tArrivalStruct,2)
            tArrivalData = [tArrivalData tArrivalStruct{j,k}];
        end
    end
    
    tArrivalMean = mean(tArrivalData);
    tArrivalMedian = median(tArrivalData);
    tArrivalSigma = std(tArrivalData);
    
    skew = skewness(tArrivalData);
    BMIN = 0;
    BMAX = 12;
    nbins = 120;
    % Uncomment to display 
    %figure; histogram(tArrivalData, nbins, 'BinLimits',[BMIN,BMAX]); title('tarrival estimate; ALL DATA');
    %xlabel(['Skewness = ' num2str(skew), ',  Mean = ' num2str(tArrivalMean), ...
    %'ms,  Median = ' num2str(tArrivalMedian), 'ms,  Sigma = ', num2str(tArrivalSigma), 'ms']);
    %set(gca, 'XLim', [BMIN BMAX])
    
end

        
