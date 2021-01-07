function [finalEst, relScore] = fuse0825(estimates, sigQualScr)
    % Given mxn array estimates where m is the number of
    % channels/estimators and n is the number of segments
    % Return a 1xn array of fused estimates
    [numChannels, numSegs] = size(estimates);
    
    finalEst = NaN(1, numSegs);
    relScore = NaN(numChannels, numSegs);
    
    consistencyScr = NaN(numChannels, numSegs);
    valScr = NaN(numChannels, numSegs);
    
    for i = 1:numSegs
        start = max([i-4, 1]);
        pastEst = estimates(:, start:i);
        
        numMatches = zeros(1, numChannels);
        for j = 1:numChannels
            numMatches(j) = sum(abs(pastEst - estimates(j, i))< 1, "all");
            consistencyScr(j, i) = max(1, numMatches(j)/numel(pastEst));
            
            if 14 < estimates(j, i) && estimates(j, i) < 20
                valScr(j, i) = 1;
            elseif 12 < estimates(j, i) && estimates(j, i) < 22
                valScr(j, i) = 0.5;
            elseif 10 < estimates(j, i) && estimates(j, i) < 24
                valScr(j, i) = 0.25;
            else
                valScr(j, i) = 0.1;
            end
        end
    end
    
    relScore = consistencyScr .* valScr .* sigQualScr;
    
    for i = 1:numSegs
        repeatVals = NaN(1, numChannels*10);
        for j = 1:numChannels
            currRepeats = round(10 * consistencyScr(j, i)^0.5 * valScr(j, i));% * (sigQualScr(j, i))^0.5);
            if isnan(currRepeats)
                continue
            end
            repeatVals(10*(j-1)+1:10*(j-1)+currRepeats) = estimates(j, i);
        end
        finalEst(i) = median(repeatVals, "omitnan");
    end
end