% %% Percent right pre fusion
% threshs = (0:1:5);
% perRight = NaN(2, numel(threshs));
% 
% allRef = [];
% allPd = [];
% allSpec = [];
% for id = 1:20
%     currRef = allEstimates(id).fuseRef;
%     currPd = allEstimates(id).pdSpecEst(4:6,:);
%     currSpec = allEstimates(id).pdSpecEst(1:3,:);
%     
%     allRef = [allRef currRef currRef currRef];
%     allPd = [allPd reshape(currPd.', [1, numel(currPd)])];
%     allSpec = [allSpec reshape(currSpec.', [1, numel(currSpec)])];
% end
% 
% totalPd = sum(~isnan(allRef) & ~isnan(allPd));
% totalSpec = sum(~isnan(allRef) & ~isnan(allSpec));
% for i = 1:numel(threshs)
%     perRight(1, i) = sum(abs(allRef - allPd) < threshs(i)) / totalPd;
%     perRight(2, i) = sum(abs(allRef - allSpec) < threshs(i)) / totalPd;
% end
% 
% perRight = round(perRight *100, 1);
% 
% coverage = NaN(1, 2);
% coverage(1) = round(totalPd / sum(~isnan(allRef)) * 100, 1);
% coverage(2) = round(totalSpec / sum(~isnan(allRef))*100, 1);
% 
% f1 = 2;
% if ishandle(f1)
%     clf(f1)
% end
% figure(f1)
% plot(threshs, perRight(1,:))
% hold on
% plot(threshs, perRight(2,:))
% xlabel("Allowed Respiratory Rate Deviation (bpm)", "FontSize", 16)
% ylabel("Percent Correct Estimates (of those that are not NaN)", "FontSize", 16)
% legend("Peak Detection Estimator ("+num2str(coverage(1))+ "% Coverage)", "Spectral Estimator ("+num2str(coverage(2))+ "% Coverage)", "FontSize", 16)
% title("Accuracy of Respiratory Rate Estimation with MI Sensors Pre Fusion", "FontSize", 20)




%% Percent right pre fusion
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
title("Accuracy of Respiratory Rate Estimation with MI Sensors Post Fusion (" + num2str(coverage) + "% Coverage)", "FontSize", 20)
