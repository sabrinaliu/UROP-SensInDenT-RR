function finalChanges = splitData(phases, segLength, fs)
% Given the phases for a patient, return a sequence of indices that
% mark the beginning of segments that are no longer than 10s and only
% include one phase.
% Also includes a index right after the end of the sequence.

    if nargin < 2
        segLength = 10; % default 10s segments
    end
    if nargin < 3
        fs = 250;
    end

    maxLength = segLength * fs;

    phaseChanges = [1 find(diff(phases))+1 numel(phases)+1]; % mark start of new phases and after end of last one

    [~, numOfChanges] = size(phaseChanges);

    finalChanges = [];
    for i = 1:numOfChanges-1
        startPhase = phaseChanges(i);
        endPhase = phaseChanges(i+1)-1;
        newChanges = startPhase:maxLength:endPhase;
        finalChanges = [finalChanges newChanges];
    end

    finalChanges = [finalChanges numel(phases)+1];
end
