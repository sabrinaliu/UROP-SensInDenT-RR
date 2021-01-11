addpath("/Users/sabrinaliu/SuperUROP/matLabFiles/heartRate/finalCode")


phaseData = struct;

for ph = 1:4
    phaseData(ph).data = [];
end

for id = 1:20
    disp(id)
    tic
    
    excp = find(patient(id).phase>0);
    phases = patient(id).phase(excp);
    
    segCuts = splitData(phases);
    numSegs = numel(segCuts)-1;
    
    
    phasePerSeg = NaN(1, numSegs);
    
    for seg = 1:numSegs
        phasePerSeg(seg) = phases(segCuts(seg));
    end
    
    currRef = allEstimates(id).refBpm;
    currRelEst = allEstimates(id).fuseBpm(1,:);
    currBayesEst = allEstimates(id).fuseBpm(2,:);
    
    for ph = 1:4
        inPhaseRef = currRef(phasePerSeg == ph);
        inPhaseRel = currRelEst(phasePerSeg == ph);
        inPhaseBayes = currBayesEst(phasePerSeg == ph);
        
        phaseData(ph).data = [phaseData(ph).data [inPhaseRef; inPhaseRel; inPhaseBayes]];
    end
end

coveragePhase = NaN(2, 4);
percentInRangePhase = NaN(2, 4);
for ph = 1:4
    currRef = phaseData(ph).data(1,:);
    currRelScr = phaseData(ph).data(2,:);
    currBayes = phaseData(ph).data(3,:);
    
    coveragePhase(1, ph) = sum(~isnan(currRelScr)) / numel(currRef);
    coveragePhase(2, ph) = sum(~isnan(currBayes)) / numel(currRef);
    
    percentInRangePhase(1, ph) = sum(abs(currRelScr - currRef) <= 5, "omitnan")/ sum(~isnan(currRelScr));
    percentInRangePhase(2, ph) = sum(abs(currBayes - currRef) <= 5, "omitnan")/ sum(~isnan(currBayes));
end


f1 = 3;
if ishandle(f1)
    clf(f1)
end
figure(f1)
for ph = 1:4
    subplot(2, 2, ph)
    blandAltman(phaseData(ph).data(1,:), phaseData(ph).data(2, :))
    title("Reliability Score Fusion for All Patients -- Phase " + num2str(ph), "FontSize", 20)
end

f1 = 4;
if ishandle(f1)
    clf(f1)
end
figure(f1)
for ph = 1:4
    subplot(2, 2, ph)
    blandAltman(phaseData(ph).data(1,:), phaseData(ph).data(3, :))
    title("Bayesian Fusion for All Patients -- Phase " + num2str(ph), "FontSize", 20)
end
