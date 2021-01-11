% %% Percent right pre fusion
% threshs = (0:1:5);
% perRight = NaN(5, numel(threshs));
% 
% allRef = [];
% allPd = [];
% allSpec = [];
% allAm = [];
% allFm = [];
% allBw = [];
% for id = 1:20
%     currRef = allEstimates(id).fuseRef;
%     currPd = allEstimates(id).pdSpecEst(4:6,:);
%     currSpec = allEstimates(id).pdSpecEst(1:3,:);
%     currAm = allEstimates(id).ppgEst(:, 1, :);
%     currFm = allEstimates(id).ppgEst(:, 2, :);
%     currBw = allEstimates(id).ppgEst(:, 3, :);
%     
%     allRef = [allRef currRef currRef currRef];
%     allPd = [allPd reshape(currPd.', [1, numel(currPd)])];
%     allSpec = [allSpec reshape(currSpec.', [1, numel(currSpec)])];
%     allAm = [allAm reshape(permute(currAm, [3, 1, 2]), [1, numel(currAm)])];
%     allFm = [allFm reshape(permute(currFm, [3, 1, 2]), [1, numel(currFm)])];
%     allBw = [allBw reshape(permute(currBw, [3, 1, 2]), [1, numel(currBw)])];
% end
% 
% totalPd = sum(~isnan(allRef) & ~isnan(allPd));
% totalSpec = sum(~isnan(allRef) & ~isnan(allSpec));
% totalAm = sum(~isnan(allRef) & ~isnan(allAm));
% totalFm = sum(~isnan(allRef) & ~isnan(allFm));
% totalBw = sum(~isnan(allRef) & ~isnan(allBw));
% for i = 1:numel(threshs)
%     perRight(1, i) = sum(abs(allRef - allPd) < threshs(i)) / totalPd;
%     perRight(2, i) = sum(abs(allRef - allSpec) < threshs(i)) / totalPd;
%     perRight(3, i) = sum(abs(allRef - allAm) < threshs(i)) / totalPd;
%     perRight(4, i) = sum(abs(allRef - allFm) < threshs(i)) / totalPd;
%     perRight(5, i) = sum(abs(allRef - allBw) < threshs(i)) / totalPd;
% end
% 
% perRight = round(perRight *100, 1);
% 
% coverage = NaN(1, 5);
% coverage(1) = round(totalPd / sum(~isnan(allRef)) * 100, 1);
% coverage(2) = round(totalSpec / sum(~isnan(allRef))*100, 1);
% coverage(3) = round(totalAm / sum(~isnan(allRef))*100, 1);
% coverage(4) = round(totalFm / sum(~isnan(allRef))*100, 1);
% coverage(5) = round(totalBw / sum(~isnan(allRef))*100, 1);
% 
% f1 = 2;
% if ishandle(f1)
%     clf(f1)
% end
% figure(f1)
% % plot(threshs, perRight(1,:))
% % hold on
% % for i = 2:5
% %     plot(threshs, perRight(i,:))
% % end
% % xlabel("Allowed Respiratory Rate Deviation (bpm)", "FontSize", 16)
% % ylabel("Percent Correct Estimates (of those that are not NaN)", "FontSize", 16)
% % legend("Peak Detection Estimator ("+num2str(coverage(1))+ "% Coverage)", "Spectral Estimator ("+num2str(coverage(2))+ "% Coverage)", "FontSize", 16)
% % title("Accuracy of Respiratory Rate Estimation with MI Sensors Pre Fusion", "FontSize", 20)
% 
% plot(threshs, perRight(3, :))
% hold on
% plot(threshs, perRight(4, :))
% plot(threshs, perRight(5, :))
% xlabel("Allowed Respiratory Rate Deviation (bpm)", "FontSize", 16)
% ylabel("Percent Correct Estimates (of those that are not NaN)", "FontSize", 16)
% legend("Amplitude Modulation Estimator ("+num2str(coverage(3))+ "% Coverage)", "Frequency Modulation Estimator ("+num2str(coverage(4))+ "% Coverage)", "Baseline Wander Estimator ("+num2str(coverage(5))+ "% Coverage)", "FontSize", 16)
% title("Accuracy of Respiratory Rate Estimation with RPPG Sensors Pre Fusion", "FontSize", 20)
% 



%% Percent right post fusion
threshs = (0:5);
perRight = NaN(1, numel(threshs));

allRef = [];
allEst = [];
for id = 1:20
    currRef = allEstimates(id).fuseRef;
    currEst = allEstimates(id).fuseEst;
    
    allRef = [allRef currRef];
    allEst = [allEst currEst];
end

totalEst = sum(~isnan(allRef) & ~isnan(allEst));
for i = 1:numel(threshs)
    perRight(i) = sum(abs(allRef - allEst) < threshs(i)) / totalEst;
end

perRight = round(perRight*100, 1);

coverage = round(totalEst / sum(~isnan(allRef)) * 100, 1);

f1 = 2;
if ishandle(f1)
    clf(f1)
end
figure(f1)
plot(threshs, perRight(1,:))
xlabel("Allowed Respiratory Rate Deviation (bpm)", "FontSize", 16)
ylabel("Percent Correct Estimates (of those that are not NaN)", "FontSize", 16)
% legend("Peak Detection Estimator ("+num2str(coverage(1))+ "% Coverage)", "Spectral Estimator ("+num2str(coverage(2))+ "% Coverage)", "FontSize", 16)
title("Accuracy of Respiratory Rate Estimation with RPPG Sensors Post Fusion (" + num2str(coverage) + "% Coverage)", "FontSize", 20)
