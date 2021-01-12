function [finalBpmEst, relScore] = spectral(inputData, filterOutStd, isPlot)
    % inputData represents 1 channel of data
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
    
%     potentialEstimates = [];
%     for i = 1:numel(winLengths)
%         if winLengths(i) > numSamples
%             continue
%         end
%         pxx = allSpectra(i).pxx;
%         bpm = allSpectra(i).bpm;
%         
%         [pkVals, pkBpms] = findpeaks(pxx, bpm, "SortStr", "descend");
%         
%         if numel(pkVals) < 1
%             continue
%         end
%         
%         if numel(potentialEstimates) == 0
%             potentialEstimates = [pkBpms pkVals ones(numel(pkBpms), 1)];
%             continue
%         end
%         
%         for j = 1:numel(pkBpms)
%             [bestDiff, bestMatchIdx] = min(abs(pkBpms(j) - potentialEstimates(:,1)));
%             
%             if bestDiff < 1
%                 potentialEstimates(bestMatchIdx, 3) = potentialEstimates(bestMatchIdx, 3) + 1;
%             else
%                 potentialEstimates = [potentialEstimates; [pkBpms(j) pkVals(j) 1]];
%             end
%         end
%         
% %         bpmEstimates(i) = pkBpms(1);
%     end
    
%     if numel(potentialEstimates) == 0
%         return
%     end
%     mostMatches = max(potentialEstimates(:,3));
%     if mostMatches < 2
%         return
%     end
%     
%     goodBpms = potentialEstimates(potentialEstimates(:,3) == mostMatches, 1:2);
%     if numel(goodBpms) == 1
%         finalBpmEst = goodBpms;
%         return
%     end
%     
%     [~, order] = sort(goodBpms(:,2));
%     goodBpms = goodBpms(order, :);
%     chosenIdx = 1;
%     for i = 2:numel(order)
%         if goodBpms(i, 1) < goodBpms(chosenIdx, 1) && goodBpms(i, 2) > .9 * goodBpms(chosenIdx, 2)
%             chosenIdx = i;
%         end
%     end
    
%     finalBpmEst = mean(bpmEstimates, "omitnan");
%     finalBpmEst = goodBpms(chosenIdx, 1);
    
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