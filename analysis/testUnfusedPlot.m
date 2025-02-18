
allRef = [];
allPdEst = [];
allPdRelScr = [];
allSpecEst = [];
allSpecRelScr = [];
for id = 1:20
    currData = allEstimates(id);
    currPd = fixDims(currData.pdSpecEst(4:6,:));
    currSpec = fixDims(currData.pdSpecEst(1:3,:));
    
    currPdRelScr = fixDims(currData.relScores(4:6,:));
    currSpecRelScr = fixDims(currData.relScores(1:3, :));
    
    allRef = [allRef currData.fuseRef currData.fuseRef currData.fuseRef];
    
    allPdEst = [allPdEst currPd];
    allSpecEst = [allSpecEst currSpec];
    
    allPdRelScr = [allPdRelScr currPdRelScr];
    allSpecRelScr = [allSpecRelScr currSpecRelScr];
end

fig = 3;
if ishandle(fig)
    clf(fig)
end
figure(fig)
ax(1) = subplot(2, 1, 1);
blandAltman(allRef, allPdEst, allPdRelScr)
title("Peak Detection Estimator for MI Sensors on All Patients")
set(gca, "FontSize", 16)
caxis([0, 1])
ax(2) = subplot(2, 1, 2);
blandAltman(allRef, allSpecEst, allSpecRelScr)
title("Spectral Approach Estimator for MI Sensors on All Patients")
set(gca, "FontSize", 16)
caxis([0, 1])
linkaxes(ax)

function output = fixDims(arr)
    output = reshape(arr.', [1, numel(arr)]);
end