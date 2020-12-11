% addpath('/Users/sabrinaliu/SuperUROP/matLabFiles/util')
% 
% allEstimates = struct;
% allCoverage = NaN(20, 3); % number for ref, est, both
% allNumSegs = NaN(20, 1);
% allNumRight = NaN(20, 2); % num fused right, num possibly right
% 
% fs = 250;
% dsFs = 5;
% 
% for id = 1:20
%     tic
%     disp(id)
%     excp = find(patient(id).phase > 0);
% 
%     time = patient(id).time(excp);
%     phase = patient(id).phase(excp);
%     resp = patient(id).ref_resp(excp);
%     coil = patient(id).coil(:, excp);
% 
%     segCuts = splitData(phase, 30);
%     numSegs = numel(segCuts) - 1;
%     
%     dsExcp = downsample(excp, fs/dsFs);
%     dsPhase = patient(id).phase(dsExcp);
%     dsResp = patient(id).ref_resp(dsExcp);
%     dsCoil = patient(id).coil(:, dsExcp);
%     
%     phaseDiffs = [0 diff(dsPhase)];
%     phaseStartIdces = [find(phaseDiffs ~= 0) numel(dsExcp)+1];
% 
%     dsSegCuts = splitData(dsPhase, 30, dsFs);
%     dsNumSegs = numel(dsSegCuts);
% 
%     refEst = NaN(2, numSegs);
%     for i = 1:numSegs
%         currRange = (segCuts(i):segCuts(i+1)-1);
%         refEst(1, i) = spectral4(resp(currRange), false);
%     end
%     refEst(2,:) = peakDetectionSliding(dsResp, phaseStartIdces, dsSegCuts);
%     fuseRefEst = fuse0825(refEst);
%     
%     pdSpecEst = NaN(6, numSegs);
%     for i = 1:numSegs
%         currRange = (segCuts(i):segCuts(i+1)-1);
%         for sensor = 1:3
%             pdSpecEst(sensor, i) = spectral4(coil(sensor, currRange));
%         end
%     end
%     for sensor = 1:3
%         pdSpecEst(sensor+3,:) = peakDetectionSliding(dsCoil(sensor,:), phaseStartIdces, dsSegCuts);
%     end
%     % fuseEst = mean(pdSpecEst, "omitnan");
%     fuseEst = fuse0825(pdSpecEst);
%     
%     allEstimates(id).numSegs = numSegs;
%     allEstimates(id).refEst = refEst;
%     allEstimates(id).fuseRef = fuseRefEst;
%     allEstimates(id).pdSpecEst = pdSpecEst;
%     allEstimates(id).fuseEst = fuseEst;
%     toc
% end

threshRight = 1;
for id = 1:20
    fuseRef = allEstimates(id).fuseRef;
    fuseEst = allEstimates(id).fuseEst;
    
    pdSpecEst = allEstimates(id).pdSpecEst;
    [numEsts, numSegs] = size(pdSpecEst);
    
    possibleRight = false(1, allEstimates(id).numSegs);
    for i = 1:numEsts
        possibleRight = possibleRight | (abs(pdSpecEst(i,:) - fuseRef) < threshRight);
    end
    
    allCoverage(id, 1) = sum(~isnan(fuseRef));
    allCoverage(id, 2) = sum(~isnan(fuseEst));
    allCoverage(id, 3) = sum(~isnan(fuseRef + fuseEst));
    
    allNumSegs(id) = numSegs;
    
    allNumRight(id, 1) = sum(abs(fuseEst - fuseRef) < threshRight);
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

% blandAltman(allRefEst, allFuseEst, allRelScr)

blandAltman(allFuseRef, allFuseEst)
title("Mean Fusion Peak Detection and Spectral Approach -- All Patients")
set(gca, "FontSize", 16)