function [finEst, relScr] = overlapPeakDetectionSliding(inputData, phaseStartIdces, segCuts, fs)
    % Given all the inputData from a channel,
    % phaseStartIdces that mark where new phases begin,
    % segCuts that mark where each estimation window begins
    % Return a finEst, 1xnumel(segCuts) array that contains the RR estimate
    % for each estimation window
    
    if nargin < 4
        fs = 5;
    end
    
    [~, numSegs] = size(segCuts);
    numSegs = numSegs-1;
    numSamples = numel(inputData);
    finEst = NaN(1, numSegs);
    relScr = NaN(1, numSegs);
    
    currPhaseStartIdx = 1;
    start = 1;
    
    unmergeEst = NaN(1, numSamples);
    while start < numSamples
        stop = start + fs*30 - 1; % 30s windows
        if stop > numSamples
            stop = numSamples;
        elseif stop >= phaseStartIdces(currPhaseStartIdx)
            stop = phaseStartIdces(currPhaseStartIdx) - 1;
            currPhaseStartIdx = currPhaseStartIdx + 1;
        end
        
        currRange = (start:stop);
        
        unmergeEst(start) = peakDetectionDouble(inputData(currRange), fs);
        
        if stop - start < 20*fs
            start = stop+1;
        else
            start = start + fs*5;
        end
    end
    
    for i = 1:numSegs
        currRange = (segCuts(1, i):segCuts(2,i));
        finEst(i) = mean(unmergeEst(currRange), "omitnan");
        relScr(i) = max(0, min(1, 1 - (std(unmergeEst(currRange), "omitnan")-1)/4));
    end
end