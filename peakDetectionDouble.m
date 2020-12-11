function [bpmEst, relScore, pkIdces] = peakDetectionDouble(unfilteredData, fs)
    bpmEst = NaN;
    relScore = NaN;
    pkIdces = [];
    
    if nargin < 2
        fs = 250;
    end
    
    if numel(unfilteredData) < 3
        return
    end
    
    filteredData = removeArtifactsRef(unfilteredData, fs);
    
%     corrData = unfilteredData - movmean(unfilteredData, fs*10);
%     [~, intervalStart, intervalEnd, relScore] = corrSliding(corrData);
%     if isnan(intervalStart)
%         intervalStart = 1;
%         intervalEnd = numel(unfilteredData);
%     end
    intervalStart = 1;
    intervalEnd = numel(unfilteredData);
    
    shortenedFilteredData = filteredData(intervalStart:intervalEnd);
    shortenedUnfilteredData = unfilteredData(intervalStart:intervalEnd);

    [bpPks, bpPkIdces] = findpeaks(shortenedFilteredData);
    [medPks, medPkIdces] = medianFilterPeaks(shortenedUnfilteredData, fs);
    
    pkMatches = nan(1, numel(bpPkIdces));
    for i = 1:numel(bpPkIdces)
        [matchDiff, matchIdx] = min(abs(bpPkIdces(i) - medPkIdces));
        if matchDiff < fs * 1.5
            pkMatches(i) = matchIdx;
        end
    end
    
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
%     relScore = assignRelScores(pkIdces, filteredData);
end

function bpmEst = assignBpmEstimate(pkIdces, fs)
    
%     bpmEst = 60 / mean(intervals, 'omitnan') * fs;

    intervals = diff(pkIdces);
    medInterval = median(intervals);
    madInterval = median(abs(intervals - medInterval));
    
    goodIntervals = intervals(abs(intervals - medInterval) < 2*madInterval);
    bpmEst = 60 / mean(goodIntervals) * fs;
end

function relScore = assignRelScores(pkIdces, inputData)
    relScore = 0;
    
    medInterval = round(median(diff(pkIdces)));
    if mod(medInterval, 2) == 0
        medInterval = medInterval + 1;
    end
    
    if numel(pkIdces) < 2
        return
    end
    
    segments = nan(numel(pkIdces), medInterval);
    halfInterval = (medInterval - 1) / 2;
    for i = 1:numel(pkIdces)
        start = pkIdces(i) - halfInterval;
        stop = pkIdces(i) + halfInterval;
        if start < 1 || stop > numel(inputData)
            continue
        end
        
        segments(i,:) = inputData(start:stop);
    end
    segments = rmmissing(segments);
    meanSeg = mean(segments);
    
    [numPks, ~] = size(segments);
    if numPks < 2
        return
    end
    
    corrs = nan(1, numPks);
    for i = 1:numPks
        currCorr = corrcoef(meanSeg, segments(i, :));
        corrs(i) = currCorr(2, 1);
    end
    
    relScore = mean(corrs);
end

function [medPks, medPkIdces] = medianFilterPeaks(unfilteredData, fs)
    numSamples = numel(unfilteredData);
    
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
    
    medFiltered = medfilt1(smoothedData, fs);
    [medPks, medPkIdces] = findpeaks(medFiltered, "MinPeakDistance", min(fs, numel(medFiltered)-2), "MinPeakProminence", 10);
end