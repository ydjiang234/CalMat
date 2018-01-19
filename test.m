close all;
clear all;
clc;
%initial the calculation class
%Use unit mm and N;
import CalMat;
A = 46759.46506;%mm2
I = 173991969.5;%mm4
L = 1750.0;%mm
Naxial = 1.0;%N
revK = 533.33E6;%Nmm/rad
ampFactor = 20.0;
backbone = load('backbone.txt');%Rotation (rad) - Moment (kNm)
%backbone(:,1) = backbone(:,1) * L;
backbone(:,2) = backbone(:,2) * 1.0e6;
targetData = load('TargetData.txt');%Displacement (mm) - Moment (kNm) 
targetData(:,1) = targetData(:,1) / L;
targetData(:,2) = targetData(:,2) * 1.0E6;

d_incr = 1.0;
cal = CalMat(A, I, L, Naxial, revK, backbone, targetData, ampFactor, d_incr);

%initial the Harmony search class
pit_range = [0, 15;
            0, 15;
            0, 15;
            0, 15;];
hms = 10;
fw_ratio = 0.01; hmcr = 0.9; par = 0.3;

HS = Harmony_Search(pit_range, hms, @cal.fit_fun, hmcr, par, fw_ratio);
max_iter = 100;
for i=1:max_iter
    HS = HS.next();
    [vectors(i,:), maxFitness(i)] = HS.Optimized();
end
%plot(maxFitness)
[temp, max_ind] =  max(maxFitness);
sprintf('With an iteration number = %i, the solutions are %0.3f %0.3f %0.3f %0.3f', max_iter, vectors(max_ind,:))