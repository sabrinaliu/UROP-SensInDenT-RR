function [finEst, relScr] = peakDetectionSliding(inputData, phaseStartIdces, segCuts, fs)
    % Calculate respiratory rate estimate for all of 30s segments in a
    % given inputData
    % Input:
    %   inputData: float vector of unfiltered data from a single channel
    %   phaseStartIdces: integer vector of indices where new
    %   phase in the dental procedure starts
    %   segCuts: integer vector of where the 30s segments start, including
    %   where the phases start
    %   fs: float storing the sampling rate
    % Output:
    %   finEst: 1xnumel(segCuts)-1 float vector of the respiratory rate
    %   estimate of each segment
    %   relScr: 1xnumel(segCuts)-1 float vector of the signal quality score
    %   corresponding to each estimate
    
    if nargin < 4
        fs = 5;
    end
    
    numSegs = numel(segCuts) - 1;
    numSamples = numel(inputData);
    finEst = NaN(1, numSegs);
    relScr = NaN(1, numSegs);
    
    currPhaseStartIdx = 1;
    start = 1;
    
    % Start Sliding window of 30s, overlap of 5s
    % In each window, apply peakDetectionDouble to get an estimate
    unmergeEst = NaN(1, numSamples);
    while start < numSamples
        stop = start + fs*30 - 1; % 30s windows
        % Prevent windows from including samples from multiple phases
        if stop > numSamples
            stop = numSamples;
        elseif stop >= phaseStartIdces(currPhaseStartIdx)
            stop = phaseStartIdces(currPhaseStartIdx) - 1;
            currPhaseStartIdx = currPhaseStartIdx + 1;
        end
        
        currRange = (start:stop);
        
        unmergeEst(start) = peakDetectionDouble(inputData(currRange), fs);
        
        % If current window <20s, then close enough to end don't need new
        % window
        if stop - start < 20*fs
            start = stop+1;
        else
            start = start + fs*5;
        end
    end
    
    % For each of the nonoverlapping 30s segments, find its estimate by
    % averaging all the sliding window estimates that overlap with this
    % segment
    for i = 1:numSegs
        currRange = (segCuts(i):segCuts(i+1)-1);
        finEst(i) = mean(unmergeEst(currRange), "omitnan");
        relScr(i) = max(0, min(1, 1 - (std(unmergeEst(currRange), "omitnan")-1)/4));
    end
end