close all;
clear all;
clc;

import Harmony_Search
pit_range = [0, 15;
            0, 15;
            0, 15;
            0, 15;];
hms = 10;
fw_ratio = 0.01; hmcr = 0.9; par = 0.3;


HS = Harmony_Search(pit_range, hms, @fit_fun, hmcr, par, fw_ratio);

max_iter = 1000;
for i=1:max_iter
    HS = HS.next();
    [vectors(i,:), maxFitness(i)] = HS.Optimized();
end
%plot(maxFitness)
[temp, max_ind] =  max(maxFitness);
sprintf('With an iteration number = %i, the solutions are %0.3f %0.3f %0.3f %0.3f', max_iter, vectors(max_ind,:))