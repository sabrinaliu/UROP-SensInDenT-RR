function [finalEst, relScore] = fuseRelScoreRR(estimates, sigQualScr)
    % Fuse estimates using reliability scores
    % Input:
    %   estimates: mxn float matrix of m respiratory rate estimates for 
    %   each of the n segments
    %   sigQualScr: mxn float matrix of corresponding signal quality scores
    %   for the estimates **currently not used for fusion
    % Output:
    %   finalEst: 1xn float matrix of final fused estimates
    %   relScores: mxn float matrix of reliability scores that were used
    
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
            % consistency score based on how many matches (estimates within
            % 1brpm of each other) from estimates from the last 4 segments
            numMatches(j) = sum(abs(pastEst - estimates(j, i))< 1, "all");
            consistencyScr(j, i) = max(1, numMatches(j)/numel(pastEst));
            
            % value score based on expectation most respiratory rates
            % within 14-20brpm and decrease in likely till get to extremes
            % of 5 and 25brpm
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
    
    % square root consistency score to increase its influence
    relScore = consistencyScr.^0.5 .* valScr;
    
    % weighted median fusion
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