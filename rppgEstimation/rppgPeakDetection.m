function [brpmEsts, rqi] = rppgPeakDetection(inputData)
    % Given a segment of one channel of RPPG data, return the respiratory
    % rate estimates based on amplitude modulation (AM), frequency
    % modulation (FM), and baseline wandering (BW)
    % Input:
    %   inputData: 1xn float vector with the original data from a segment
    %   of RPPG sensor
    % Output:
    %   brpmEsts: 1x3 vector with the AM, FM, BW estimates in that order
    %   rqi: 1x3 vector with the respiratory quality indices for each est
    
    brpmEsts = [];
    rqi = [];

    fs = 250;
    
    % Smooth and filter data, identify peaks
    semiCleanData = cleaning(inputData); 
    [semiCleanPks, semiCleanPkIdces] = findpeaks(semiCleanData);
    
    % For each nonoverlapping, 10s segment, filter and identify heart rate
    % peaks using peak detection
    fullCleanPkIdces = [];
    currIdx = 1;
    allRREst = [];
    allRRPdRelScore = [];
    while currIdx <= numel(inputData)
        if numel(inputData) - currIdx < 20 * fs
            % If less than 20s left, just include the entire rest of the
            % sequence as part of this window's estimate
            % Prevents very short windows that are <10s at the end
            stop = numel(inputData);
        else
            stop = currIdx + 10*fs - 1;
        end
        
        currCleanData = removeArtifactsHR(inputData(currIdx:stop)); % Filtering from HR estimation
        [currRREst, currRRPdRelScore, currPkIdces] = peakDetectionHR(currCleanData); % peak detection from HR estimation
        fullCleanPkIdces = [fullCleanPkIdces currPkIdces+currIdx-1];
        allRREst = [allRREst currRREst];
        allRRPdRelScore = [allRRPdRelScore currRRPdRelScore];
        currIdx = stop + 1; % start next segment after previous one ended
    end
    
    % Compare peaks in the two filtered versions
    % If peaks are .2s away, they are considered the same peak
    actualSemiCleanPkIdces = [];
    for i = 1:numel(fullCleanPkIdces)
        potentialMatches = abs(semiCleanPkIdces - fullCleanPkIdces(i)) < .2*fs;
        potentialPks = semiCleanPks(potentialMatches);
        potentialPkIdces = semiCleanPkIdces(potentialMatches);

        [~, bestPkIdxOfPotential] = max(potentialPks);
        bestPkIdx = potentialPkIdces(bestPkIdxOfPotential);

        actualSemiCleanPkIdces = [actualSemiCleanPkIdces bestPkIdx];
    end
    
    % Only keep the unique peaks
    actualSemiCleanPkIdces = unique(actualSemiCleanPkIdces);
    if numel(actualSemiCleanPkIdces) < 2
        return
    end
    
    % Locate the valley points in between the peaks
    actualSemiCleanVlIdces = [];
    for i = 1:numel(actualSemiCleanPkIdces)
        if i == 1
            start = 1;
        else
            start = actualSemiCleanPkIdces(i-1);
        end

        stop = actualSemiCleanPkIdces(i);

        [~, bestVlIdxOfRange] = min(semiCleanData(start:stop));
        actualSemiCleanVlIdces = [actualSemiCleanVlIdces bestVlIdxOfRange + start - 1];
    end
    
    if numel(actualSemiCleanPkIdces) < 3
        return
    end
    
    % Based on peaks and valleys, identify amplitude, frequency, and
    % baseline values for each cardiac cycle
    actualSemiCleanPkTime = actualSemiCleanPkIdces / fs;
    am = [actualSemiCleanPkTime; semiCleanData(actualSemiCleanPkIdces) - semiCleanData(actualSemiCleanVlIdces)];
    fm = [actualSemiCleanPkTime(1:end-1); actualSemiCleanPkIdces(2:end) - actualSemiCleanPkIdces(1:end-1)];
    bw = [actualSemiCleanPkTime; semiCleanData(actualSemiCleanPkIdces)];
    
    % Clean the modulations, removing anything with sharp jumps
    am = cleanMods(am);
    fm = cleanMods(fm);
    bw = cleanMods(bw);
    
    % Interpolate the values for the modulations so that they are evenly
    % spaced and sampled at 4Hz
    newFs = 4;
    indexRange = (min(actualSemiCleanPkTime):1/newFs: max(actualSemiCleanPkTime));
        
    amMid = spline(am(1,:), am(2,:), indexRange);
    fmMid = spline(fm(1,:), fm(2,:), indexRange);
    bwMid = spline(bw(1,:), bw(2,:), indexRange);
    
    % Estimate the respiratory rate using correlation
    [brpmEsts, rqi] = bpmEstsCorr(amMid, fmMid, bwMid, newFs);
end

function [bpmEsts, rqi] = bpmEstsCorr(am, fm, bw, newFs)
    % Given evenly sampled modulations at the sampling frequency newFs,
    % return the corresponding respiratory rates in a 1x3 float vector
    % bpmEsts and the corresponding respiratory quality indices in a 1x3
    % float vector rqi
    
    bpmEsts = NaN(1, 3);
    rqi = NaN(1, 3);
    
    % Allowable respiratory rates
    minBpm = 5;
    maxBpm = 25;
    
    [bpmEsts(1), rqi(1)] = getCorrAndRqi(am);
    [bpmEsts(2), rqi(2)] = getCorrAndRqi(fm);
    [bpmEsts(3), rqi(3)] = getCorrAndRqi(bw);
    
    function [currBpmEst, currRqi] = getCorrAndRqi(modulation)
        % Using autocorrelation, find an estimate for the respiratory rate
        
        currBpmEst = NaN;
        currRqi = NaN;
        
        N = numel(modulation); % number of samples
    
        % Given min and max brpm, find the range of possible number of
        % samples that can be in a respiratory cycle
        rrRange = (floor(60*newFs/maxBpm) : ceil(60*newFs/minBpm));
        
        modulation = modulation - mean(modulation);

        % Find the autocorrelation values
        corrVals = NaN(1, numel(rrRange));
        for k = rrRange
            if k > N
                break
            end
            laggedMod = NaN(1, numel(modulation));
            laggedMod(k+1:end) = modulation(1:N-k);

            corrVals(k - rrRange(1) + 1) = sum(modulation .* laggedMod, "omitnan") / sum(modulation.^2);
        end
        
        % Identify the best shift by finding peaks in the autocorrelation
        [currRqis, bestLags] = findpeaks(corrVals, rrRange);
        if ~isnan(bestLags)
            bestBpms = 60 * newFs ./ bestLags;
            [chosenBpm, ~] = max(bestBpms);
            
            % Due to peaks of larger magnitude, may estimate multiples of
            % the actual brpm, include those when considering reliability
            % score (rqi)
            relevantBpms = bestBpms;
            for i = 1:numel(bestBpms)
                ratio = chosenBpm / bestBpms(i);
                if abs(ratio - round(ratio)) > 0.1
                    relevantBpms(i) = NaN;
                end
            end
            currRqis(isnan(relevantBpms)) = [];
            
            currRqi = sum(currRqis);
            currBpmEst = chosenBpm;
        end
    end
    
end

function newMod = cleanMods(mod)
    % Clean given modulations mod by removing any sharp jumps/ drops
    
    newMod = mod;
    
    diffMod = diff(mod, 1, 2);
    derivs = diffMod(2,:) ./ diffMod(1,:);
    rmsDeriv = rms(derivs);
    
    badDerivIdces = find(abs(derivs) > 1.5*rmsDeriv);
    newMod(:, badDerivIdces) = NaN;
    newMod(:, badDerivIdces+1) = NaN;
    
    newMod = rmmissing(newMod, 2);
end

function finalData = cleaning(inputData)
    % Smooth and clean the raw RPPG data for respiratory rate estimation
    
    fs = 250;
    numSamples = numel(inputData);

    miniWin = 5 * fs; % 1.5s window length
    smoothedData = nan(1, numSamples);
    for winStart = 1:miniWin:numSamples
        winEnd = min(winStart+miniWin-1, numSamples);
        currSection = inputData(winStart:winEnd);

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
        smoothedData(winStart:winEnd) = newSection;
    end

    smoothedData = sgolayfilt(smoothedData, 2, 45); % smoothing
    
    smoothedData = smoothedData - mean(smoothedData);
    [b, a] = butter(5, 3/fs*2); % low pass filter to be within respiratory rate range
    finalData = filtfilt(b, a, smoothedData);
end
