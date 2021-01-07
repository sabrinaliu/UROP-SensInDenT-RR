% addpath('/Users/sabrinaliu/SuperUROP/matLabFiles/util')
% 
% allEstimatesOverlap = struct;
% allCoverageOverlap = NaN(20, 3); % number for ref, est, both
% allNumSegsOverlap = NaN(20, 1);
% allNumRightOverlap = NaN(20, 2); % num fused right, num possibly right

fs = 250;
dsFs = 5;
idRange = 10:10;

for id = idRange
    tic
    disp(id)
    excp = find(patient(id).phase > 0);

    time = patient(id).time(excp);
    phase = patient(id).phase(excp);
    resp = patient(id).ref_resp(excp);
    coil = patient(id).coil(:, excp);

    segCuts = splitDataOverlap(phase, 30, 15);
    [~, numSegs] = size(segCuts);
    numSegs = numSegs - 1;
    
    dsExcp = downsample(excp, fs/dsFs);
    dsPhase = patient(id).phase(dsExcp);
    dsResp = patient(id).ref_resp(dsExcp);
    dsCoil = patient(id).coil(:, dsExcp);
    
    phaseDiffs = [0 diff(dsPhase)];
    phaseStartIdces = [find(phaseDiffs ~= 0) numel(dsExcp)+1];

    dsSegCuts = splitDataOverlap(dsPhase, 30, 15, dsFs);

    refEst = NaN(2, numSegs);
    refRelScr = NaN(2, numSegs);
    for i = 1:numSegs
        currRange = (segCuts(1, i):segCuts(2, i));
        [refEst(1, i), refRelScr(1, i)] = spectral4(resp(currRange), false);
    end
     [currRefEst, currRelScr] = overlapPeakDetectionSliding(dsResp, phaseStartIdces, dsSegCuts);
     refEst(2,1:numel(currRefEst)) = currRefEst; 
     refRelScr(2,1:numel(currRefEst)) = currRelScr;
    fuseRefEst = fuse0825(refEst, refRelScr);
    
    pdSpecEst = NaN(6, numSegs);
    pdSpecRelScr = NaN(6, numSegs);
    for i = 1:numSegs
        currRange = (segCuts(i):segCuts(2, i));
        for sensor = 1:3
            [pdSpecEst(sensor, i), pdSpecRelScr(sensor, i)] = spectral4(coil(sensor, currRange));
        end
    end
    for sensor = 1:3
         [currPdEst, currPdRelScr] = overlapPeakDetectionSliding(dsCoil(sensor,:), phaseStartIdces, dsSegCuts);
         pdSpecEst(sensor+3,1:numel(currPdEst)) = currPdEst;
         pdSpecRelScr(sensor+3,1:numel(currPdEst)) = currPdRelScr;
    end
    fuseEst = fuse0825(pdSpecEst, pdSpecRelScr);
    
    allEstimatesOverlap(id).numSegs = numSegs;
    allEstimatesOverlap(id).refEst = refEst;
    allEstimatesOverlap(id).refRelScr = refRelScr;
    allEstimatesOverlap(id).fuseRef = fuseRefEst;
    allEstimatesOverlap(id).pdSpecEst = pdSpecEst;
    allEstimatesOverlap(id).pdSpecRelScr = pdSpecRelScr;
    allEstimatesOverlap(id).fuseEst = fuseEst;
    toc
end

threshRight = 1;
for id = idRange
    fuseRef = allEstimatesOverlap(id).fuseRef;
    fuseEst = allEstimatesOverlap(id).fuseEst;
    
    pdSpecEst = allEstimatesOverlap(id).pdSpecEst;
    [numEsts, numSegs] = size(pdSpecEst);
    
    possibleRight = false(1, allEstimatesOverlap(id).numSegs);
    for i = 1:numEsts
        possibleRight = possibleRight | (abs(pdSpecEst(i,:) - fuseRef) < threshRight);
    end
    
    allCoverageOverlap(id, 1) = sum(~isnan(fuseRef));
    allCoverageOverlap(id, 2) = sum(~isnan(fuseEst));
    allCoverageOverlap(id, 3) = sum(~isnan(fuseRef + fuseEst));
    
    allNumSegsOverlap(id) = numSegs;
    
    allNumRightOverlap(id, 1) = sum(abs(fuseEst - fuseRef) < threshRight);
    allNumRightOverlap(id, 2) = sum(possibleRight);
end

allFuseRef = [];
allFuseEst = [];
for id = idRange
    allFuseRef = [allFuseRef allEstimatesOverlap(id).fuseRef];
    allFuseEst = [allFuseEst allEstimatesOverlap(id).fuseEst];
end

fig2 = 4;
if ishandle(fig2)
    clf(fig2)
end
figure(fig2)
timeStart = segCuts(1,1:end-1);
timeStart = time(timeStart);
stairs(timeStart, allFuseRef)
hold on
stairs(timeStart, allFuseEst)
xlabel("Time (s)")
ylabel("Amplitude (a.u.)")
title("Patient 10 with Overlap")
set(gca, "FontSize", 16)

fig = 1;
if ishandle(fig)
    clf(fig)
end
figure(fig)

% blandAltman(allRefEst, allFuseEst, allRelScr)

blandAltman(allFuseRef, allFuseEst)
title("Fusion Peak Detection Overlap -- Patient 10")
set(gca, "FontSize", 16)