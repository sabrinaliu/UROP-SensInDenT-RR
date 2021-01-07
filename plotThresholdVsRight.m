
threshs = (0:0.5:5);
perRight = NaN(1, numel(threshs));

for idx = 1:numel(threshs)
    currThresh = threshs(idx);
    
    allTotal = 0;
    allRight = 0;
    for id = 1:20
        currRef = allEstimates(id).fuseRef;
        currEst = allEstimates(id).fuseEst;
        
        total = sum(~isnan(currRef) & ~isnan(currEst));
        right = sum(abs(currRef - currEst) < currThresh);
        
        allTotal = allTotal + total;
        allRight = allRight + right;
    end
    perRight(idx) = allRight / allTotal;
end

coverage = round(sum(~isnan(currRef) & ~isnan(currEst)) / sum(~isnan(currRef)) * 100, 1);

f = 4;
if ishandle(f)
    clf(f)
end
figure(f)
plot(threshs, perRight)
xlabel("Allowed Respiratory Rate Deviation (brpm)")
ylabel("Percent Correst Estimates")
title("Accuracy of Respiratory Rate Estimation with MI Sensors (Coverage " + num2str(coverage) + "%)")

set(gca, "FontSize", 20)

