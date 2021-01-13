function finalData = removeArtifacts(inputData, fs)
% Given the original signal (i.e. patient(id).coil(1,:)), return a new signal
% that has been  filtered for artifacts for respiratory rate estimation
% Input:
%   inputData: 1xn float vector with the original signal
%   fs: float with the sampling frequency (default 250Hz)
% Output:
%   filtered: 1xn float vector with filtered signal
%   numArtifacts: integer with number of points identified as a sharp spike
%   artifact
%   hasArtifact: 1xn boolean vector that identifies if each point was part
%   of an artifact that needed to be removed and interpolated

    if nargin < 2
        fs = 250;
    end
    fmin = 5/60;
    fmax = 25/60;
    Nmax = ceil(fs/fmin);
    
    numSamples = numel(inputData);
    
    % Remove all points that correspond to sharp spikes and interpolate
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
    
    % Bandpass filter to 5-25brpm
    filtered = smoothedData;
    if numel(filtered) > 30
        Wn = [fmin, fmax] ./ fs .* 2;
        [b, a] = butter(3, Wn);
        filtered = filtfilt(b, a, filtered);
    end
    
    % Subtract mean to center around 0 and remove baseline wandering
    baseline = movmean(filtered, floor(1.5*Nmax), 'omitnan');
    finalData = filtered - baseline;
    
end