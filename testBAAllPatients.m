addpath('/Users/sabrinaliu/SuperUROP/matLabFiles/respiratoryRate/finalCode/util')

allEstimates = struct;
allCoverage = NaN(20, 3); % number for ref, est, both
allNumSegs = NaN(20, 1);
allNumRight = NaN(20, 2); % num fused right, num possibly right

fs = 250;
dsFs = 5;

for id = 1:20
    tic
    disp(id)
    excp = find(patient(id).phase > 0);

    time = patient(id).time(excp);
    phase = patient(id).phase(excp);
    resp = patient(id).ref_resp(excp);
    coil = patient(id).coil(:, excp);

    segCuts = splitData(phase, 30);
    numSegs = numel(segCuts) - 1;
    
    dsExcp = downsample(excp, fs/dsFs);
    dsPhase = patient(id).phase(dsExcp);
    dsResp = patient(id).ref_resp(dsExcp);
    dsCoil = patient(id).coil(:, dsExcp);
    
    phaseDiffs = [0 diff(dsPhase)];
    phaseStartIdces = [find(phaseDiffs ~= 0) numel(dsExcp)+1];

    dsSegCuts = splitData(dsPhase, 30, dsFs);
    dsNumSegs = numel(dsSegCuts);

    refEst = NaN(2, numSegs);
    refRelScr = NaN(2, numSegs);
    for i = 1:numSegs
        currRange = (segCuts(i):segCuts(i+1)-1);
        [refEst(1, i), refRelScr(1, i)] = spectral(resp(currRange), false);
    end
    [refEst(2,:), refRelScr(2,:)] = peakDetectionSliding(dsResp, phaseStartIdces, dsSegCuts);
    fuseRefEst = fuseRelScoreRR(refEst, refRelScr);
    
    pdSpecEst = NaN(6, numSegs);
    pdSpecRelScr = NaN(6, numSegs);
    for i = 1:numSegs
        currRange = (segCuts(i):segCuts(i+1)-1);
        for sensor = 1:3
            [pdSpecEst(sensor, i), pdSpecRelScr(sensor, i)] = spectral(coil(sensor, currRange));
        end
    end
    for sensor = 1:3
        [pdSpecEst(sensor+3,:), pdSpecRelScr(sensor+3,:)] = peakDetectionSliding(dsCoil(sensor,:), phaseStartIdces, dsSegCuts);
    end
    [fuseEst, relScores] = fuseRelScoreRR(pdSpecEst, pdSpecRelScr);
    
    allEstimates(id).numSegs = numSegs;
    allEstimates(id).refEst = refEst;
    allEstimates(id).refRelScr = refRelScr;
    allEstimates(id).fuseRef = fuseRefEst;
    allEstimates(id).pdSpecEst = pdSpecEst;
    allEstimates(id).pdSpecRelScr = pdSpecRelScr;
    allEstimates(id).fuseEst = fuseEst;
    allEstimates(id).relScores = relScores;
    toc
end

for id = 1:20
    fuseRef = allEstimates(id).fuseRef;
    fuseEst = allEstimates(id).fuseEst;
    
    pdSpecEst = allEstimates(id).pdSpecEst;
    [numEsts, numSegs] = size(pdSpecEst);
    
    possibleRight = false(1, allEstimates(id).numSegs);
    for i = 1:numEsts
        possibleRight = possibleRight | (abs(pdSpecEst(i,:) - fuseRef) < 1);
    end
    
    allCoverage(id, 1) = sum(~isnan(fuseRef));
    allCoverage(id, 2) = sum(~isnan(fuseEst));
    allCoverage(id, 3) = sum(~isnan(fuseRef + fuseEst));
    
    allNumSegs(id) = numSegs;
    
    allNumRight(id, 1) = sum(abs(fuseEst - fuseRef) < 1);
    allNumRight(id, 2) = sum(possibleRight);
end

allFuseRef = [];
allFuseEst = [];
for id = 1:20
    allFuseRef = [allFuseRef allEstimates(id).fuseRef];
    allFuseEst = [allFuseEst allEstimates(id).fuseEst];
end

fig = 1;
if ishandle(fig)
    clf(fig)
end
figure(fig)

blandAltman(allFuseRef, allFuseEst)
title("Respiratory Rate Post Fusion Estimation with MI Sensors")
set(gca, "FontSize", 20)
