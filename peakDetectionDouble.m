function [bpmEst, relScore, pkIdces] = peakDetectionDouble(unfilteredData, fs)
    % Subroutine used by peakDetectionSliding
    % Applies peak detection to signal filtered 2 different ways (using
    % bandpass filter and using median filter) and compares the two to
    % remove extraneous peaks
    % Input:
    %   unfilteredData: 1xn float vector of unfiltered data from a window
    %   fs: float for the sampling rate
    % Output:
    %   bpmEst: float for the respiratory rate estimate for this window
    %   relScore: float for the signal quality score for that estimate
    %   pkIdces: float vector with indices of the peaks identified
    
    bpmEst = NaN;
    relScore = NaN;
    pkIdces = [];
    
    if nargin < 2
        fs = 250;
    end
    
    if numel(unfilteredData) < 3
        return
    end
    
    filteredData = removeArtifacts(unfilteredData, fs); % normal filtering with bandmass

    [bpPks, bpPkIdces] = findpeaks(filteredData);
    [medPks, medPkIdces] = medianFilterPeaks(unfilteredData, fs); % filtering with median
    
    % Use median peaks to cross check if the bandpass peaks are extraneous
    % or not
    pkMatches = nan(1, numel(bpPkIdces));
    for i = 1:numel(bpPkIdces)
        [matchDiff, matchIdx] = min(abs(bpPkIdces(i) - medPkIdces));
        if matchDiff < fs * 1.5 % consider a match if within 1.5s of each other
            pkMatches(i) = matchIdx;
        end
    end
    
    % if multiple bandpass peaks got matched to the same median peak, then
    % only match the median peak to the tallest (greatest magnitude)
    % bandpass peak
    for i = 1:numel(bpPkIdces)
        if isnan(pkMatches(i))
            continue
        end
        
        sameMatchIdces = find(pkMatches == pkMatches(i));
        if numel(sameMatchIdces) == 1
            continue
        end
        
        sameMatchPkHeights = bpPks(sameMatchIdces);
        [~, maxMatchIdxOfSameMatch] = max(sameMatchPkHeights);
        
        maxMatchIdx = sameMatchIdces(maxMatchIdxOfSameMatch);
        
        for j = 1:numel(sameMatchIdces)
            if sameMatchIdces(j) ~= maxMatchIdx
                pkMatches(sameMatchIdces(j)) = NaN;
            end
        end
    end
    
    pkIdces = bpPkIdces(~isnan(pkMatches));
    bpmEst = assignBpmEstimate(pkIdces, fs);
end

function bpmEst = assignBpmEstimate(pkIdces, fs)
    % Get respiratory rate estimate using a MAD process, if there are peaks
    % that are too far or too close together, ignore them
    intervals = diff(pkIdces);
    medInterval = median(intervals);
    madInterval = median(abs(intervals - medInterval));
    
    goodIntervals = intervals(abs(intervals - medInterval) < 2*madInterval);
    bpmEst = 60 / mean(goodIntervals) * fs;
end

function [medPks, medPkIdces] = medianFilterPeaks(unfilteredData, fs)
    % Apply artifact removal process with median filter instead of bandpass
    % filter
    % Afterwards, identifies peaks: their magnitudes (medPks) and their
    % indices (medPkIdces)
    
    numSamples = numel(unfilteredData);
    
    % Identify sharp spikes in small 10s windows and interpolate
    miniWin = 10 * fs; % 10s window length
    for winStart = 1:miniWin:numSamples
        winEnd = min(winStart+miniWin-1, numSamples);
        if numSamples - winEnd < 5
            winEnd = numSamples;
        end
        currSection = unfilteredData(winStart:winEnd);

        currDeriv = diff(currSection);
        currRms = rms(currDeriv);

        removeIdcs = currDeriv < -1.5*currRms | 1.5*currRms < currDeriv;
        currSection(removeIdcs) = NaN;
        currSection(removeIdcs+1) = NaN;
        currSampPts = (winStart:winEnd);
        currSampPts(isnan(currSection)) = NaN;

        currSection = rmmissing(currSection);
        currSampPts = rmmissing(currSampPts);
        
        newSection = pchip(currSampPts, currSection, (winStart:winEnd)); % interpolate pts
        newRef(winStart:winEnd) = newSection;
        
        if winEnd == numSamples
            break
        end
    end
    
    smoothedData = newRef;
    if fs > 100
        smoothedData = sgolayfilt(newRef, 2, 45); % smoothing
    end
    
    % Median filter to help further smooth the data and remove high
    % frequency noise
    medFiltered = medfilt1(smoothedData, fs);
    
    % Only keep peaks that are at least 1s apart from each other (expect
    % respiratory cycles to be a certain length) and have a peak prominence
    % of at least 10 (avoid small peaks from noise)
    [medPks, medPkIdces] = findpeaks(medFiltered, "MinPeakDistance", min(fs, numel(medFiltered)-2), "MinPeakProminence", 10);
end