function refBpm = getRefBpm(inputData, segCuts)
% Given the ecg signal as inputData and indices segCuts for segment cuts,
% return an array that contains the average frequency in each of the
% segments to be the reference

    fs = 250;
    [~, rIdcs, ~] = pan_tompkin(inputData, fs, 0);

    refBpm = zeros(1, numel(segCuts)-1);
    for i = 1:numel(segCuts) - 1
        currRIdcs = rIdcs(segCuts(i) <= rIdcs & rIdcs < segCuts(i+1));
        refBpm(i) =  fs / mean(diff(currRIdcs)) * 60;

        if isnan(refBpm(i)) & ~isnan(refBpm(i-1))
            refBpm(i) = refBpm(i-1);
        elseif isnan(refBpm(i))
            refBpm(i) = 80;
        end
    end

end
