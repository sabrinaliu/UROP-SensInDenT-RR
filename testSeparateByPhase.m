
for id = 1:20

    excp = find(patient(id).phase > 0);

    time = patient(id).time(excp);
    phase = patient(id).phase(excp);
    resp = patient(id).ref_resp(excp);
    coil = patient(id).coil(:, excp);

    segCuts = splitData(phase, 30);
    numSegs = numel(segCuts) - 1;
    
    allEstimates(id).phaseSegs = phase(segCuts(1:end-1));
end

byPhase = struct;
for phase = 1:4
    byPhase(phase).ref = [];
    byPhase(phase).est = [];
    for id = 1:20
        currPhaseSegs = allEstimates(id).phaseSegs;
        newRef = allEstimates(id).fuseRef(currPhaseSegs == phase);
        newEst = allEstimates(id).fuseEst(currPhaseSegs == phase);
        
        byPhase(phase).ref = [byPhase(phase).ref newRef];
        byPhase(phase).est = [byPhase(phase).est newEst];
    end
end

perCoverage = NaN(1, 4);
perNumRight = NaN(1, 4);
threshRight = 1;
for phase = 1:4
    currRef = byPhase(phase).ref;
    currEst = byPhase(phase).est;
    total = sum(~isnan(currRef));
    perCoverage(phase) = sum(~isnan(currEst) & ~isnan(currRef)) / total;
    perNumRight(phase) = sum(abs(currEst - currRef) < threshRight) / total;
end

fig = 2;
if ishandle(fig)
    clf(fig)
end
figure(fig)
for phase = 1:4
    ax(phase) = subplot(2, 2, phase);
    blandAltman(byPhase(phase).ref, byPhase(phase).est)
    title("Fusion Peak Detection and Spectral Approach All Patients Phase " + num2str(phase))
    set(gca, "FontSize", 16)
end
linkaxes(ax)