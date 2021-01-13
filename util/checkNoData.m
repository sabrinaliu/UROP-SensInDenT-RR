function isNoData = checkNoData(inputData)
% Given a raw signal from multiple sensors, check if it has any data in any
% of the sensor channels

    [~, numPts] = size(inputData);

    if numPts <= 1
        isNoData = true;
        return
    end

    isNoData = sum(diff(inputData, 1, 2), 'all') == 0;
end
