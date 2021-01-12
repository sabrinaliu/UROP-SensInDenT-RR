function finalData = removeArtifacts(inputData, fs)
    % Removes sharp spikes, smooths, filters the data from one channel
    
    if nargin < 2
        fs = 250;
    end
    fmin = 5/60;
    fmax = 25/60;
    Nmax = ceil(fs/fmin);
    
    numSamples = numel(inputData);
    
    miniWin = 10 * fs; % 10s window length
    for winStart = 1:miniWin:numSamples
        winEnd = min(winStart+miniWin-1, numSamples);
        if numSamples - winEnd < 5
            winEnd = numSamples;
        end
        
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
        
        if numel(currSection) < 3
            inputData(winStart:winEnd) = inputData(winStart:winEnd);
            continue
        end
        newSection = pchip(currSampPts, currSection, (winStart:winEnd)); % interpolate pts
        inputData(winStart:winEnd) = newSection;
        
        if winEnd == numSamples
            break
        end
    end
    
    smoothedData = inputData;
    if fs > 100
        fl = 21;
        smoothedData = sgolayfilt(inputData, 2, fl); % smoothing
    end
    
    filtered = smoothedData;
    if numel(filtered) > 30
        Wn = [fmin, fmax] ./ fs .* 2;
        [b, a] = butter(3, Wn);
        filtered = filtfilt(b, a, filtered);
    end
    
    baseline = movmean(filtered, floor(1.5*Nmax), 'omitnan');
    finalData = filtered - baseline;
    
end