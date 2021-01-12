allEstimates = struct;
for id = 1:20
    currData = allEstimatesboth(id);
    currRef = currData.fuseRef;
    currPpg = currData.ppgEst;
    currRqi = currData.ppgRqi;
    
    numSamps = numel(currRef);
    
    currAm = permute(currPpg(:, 1, :), [1, 3, 2]);
    currFm = permute(currPpg(:, 2, :), [1, 3, 2]);
    currBw = permute(currPpg(:, 3, :), [1, 3, 2]);
    
    [finalEst, relScore] = fuseRelScoreRR([currAm; currFm; currBw], ones(9, numSamps));
    
    allEstimates(id).amFmBw = [currAm; currFm; currBw];
    allEstimates(id).fuseEst = finalEst;
    allEstimates(id).fuseRef = currRef;
    allEstimates(id).relScore = relScore;
end

allRef = [];
allEst = [];
for id = 1:20
    allRef = [allRef allEstimates(id).fuseRef];
    allEst = [allEst allEstimates(id).fuseEst];
end

fig = 1;
if ishandle(fig)
    clf(fig)
end
figure(fig)

blandAltman(allFuseRef, allFuseEst)

% blandAltman(allFuseRef, allFuseEst, allThird)
title("Respiratory Rate Post Fusion Estimation with RPPG")
set(gca, "FontSize", 20)