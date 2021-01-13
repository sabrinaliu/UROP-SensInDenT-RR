function [finalBpmEst, relScore] = spectral(inputData, filterOutStd, isPlot)
    % Given data from one segment, use spectral approach to estimate
    % respiratory rate
    % Input:
    %   inputData: 1xn float vector that includes data from one segment
    %   filterOutStd: boolean whether should omit estimates with high
    %   variation
    %   isPlot: boolean if should plot results
    % Output:
    %   finalBpmEst: float that is the estimate of current respiratory rate
    %   in brpm
    %   relScore: float that is the reliability score of the returned
    %   estimate
    
    finalBpmEst = NaN;
    relScore = NaN;
    
    if nargin < 2
        filterOutStd = true;
        isPlot = false;
    end
    if nargin < 3
        isPlot = false;
    end
    
    if sum(isnan(inputData)) == numel(inputData)
        finalBpmEst = NaN;
        relScore = 0;
        return
    end
    
    fs = 250;
    [~, numSamples] = size(inputData);
    if numSamples > 1
        inputData = inputData';
        [numSamples, ~] = size(inputData);
    end
    
    winLengths = [30, 20, 15, 12] .* fs;
    
    allSpectra = struct;
    for i = 1:numel(winLengths)
        if winLengths(i) > numSamples
            continue
        end
        
        [pxx, bpm] = spectra(inputData, winLengths(i));
        
        allSpectra(i).pxx = pxx;
        allSpectra(i).bpm = bpm;
    end
    
    bpmEstimates = NaN(1, numel(winLengths));
    for i = 1:numel(winLengths)
        if winLengths(i) > numSamples
            continue
        end
        pxx = allSpectra(i).pxx;
        bpm = allSpectra(i).bpm;
        
        [~, pkBpms] = findpeaks(pxx, bpm, "SortStr", "descend", "NPeaks", 1);
        if numel(pkBpms) == 1
            bpmEstimates(i) = pkBpms;
        end
    end
    
    if ~filterOutStd || std(bpmEstimates, "omitnan") < 4
        relScore = max(0, min(1, 1 - (std(bpmEstimates, "omitnan")-1)/3));
        finalBpmEst = mean(bpmEstimates, "omitnan");
    end
    
    if isPlot
        fig = 5;
        if ishandle(fig)
            clf(fig)
        end
        figure(fig)
        ax = NaN;
        for i = 1:numel(winLengths)
            pxx = allSpectra(i).pxx;
            bpm = allSpectra(i).bpm;

            ax(i) = subplot(numel(winLengths), 1, i);
            plot(bpm, pxx)
            hold on
            if ~isnan(bpmEstimates(i))
                scatter(bpmEstimates(i), 0)
            end
        end
    end
end

function [pxx, bpm] = spectra(data, winLength)
    fs = 250;
    
    % Compute spectra
    data = data(end-winLength+1:end, :);
    data = data - movmean(data, fs*3);
    
    [numSamples, ~] = size(data);
    
    [pxx, f] = pwelch(data, numSamples, [], [], fs);
    pxx = sum(pxx,2);
    pxx = 10 .* log10(pxx);
    
    bpm = f * 60;
    pxx = pxx(5 < bpm & bpm < 25);
    bpm = bpm(5 < bpm & bpm < 25);
end