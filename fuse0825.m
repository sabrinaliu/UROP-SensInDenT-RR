function [finalEst, relScore] = fuse0825(estimates)
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
            
            if estimates(j, i) < 10 || estimates(j, i) > 20
                valScr(j, i) = .5;
            else
                valScr(j, i) = 1;
            end
        end
    end
    
    relScore = consistencyScr .* valScr;
    
    for i = 1:numSegs
        repeatVals = NaN(1, numChannels*10);
        for j = 1:numChannels
            currRepeats = round(10 * consistencyScr(j, i) * valScr(j, i));
            repeatVals(10*(j-1)+1:10*(j-1)+currRepeats) = estimates(j, i);
        end
        finalEst(i) = median(repeatVals, "omitnan");
    end
end