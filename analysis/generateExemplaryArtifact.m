%% Exemplary artifact MI
id = 20;
seg = 10;
sensor = 3;

excp = find(patient(id).phase>0);
time = patient(id).time(excp);
ecg = patient(id).ref_ecg(excp);
coil = patient(id).coil(:, excp);
phases = patient(id).phase(excp);

segCuts = splitData(phases, 30);

currRange = (segCuts(seg):segCuts(seg+1)-1);

filtered = removeArtifactsRef(coil(sensor, currRange));

f2 = 4;
if ishandle(f2)
    clf(f2)
end
figure(f2)
subplot(2, 1, 1)
plot(time(currRange), coil(sensor, currRange))
xlabel("Time (s)")
ylabel("Amplitude (a.u.)")
title("Patient " + num2str(id) + " Sensor " + num2str(sensor) + " Segment " + num2str(seg) + " Original")
set(gca, "FontSize", 20)

subplot(2, 1, 2)
plot(time(currRange), filtered)
xlabel("Time (s)")
ylabel("Amplitude (a.u.)")
title("Patient " + num2str(id) + " Sensor " + num2str(sensor) + " Segment " + num2str(seg) + " Filtered")
set(gca, "FontSize", 20)