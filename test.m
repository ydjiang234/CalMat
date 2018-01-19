close all;
clear all;
clc;
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
cal = CalMat(A, I, L, Naxial, revK, backbone, targetData, ampFactor);

%filePath = 'ModifiedCPTemplate.tcl';
[output, energy, fitness] = cal.Analyze(10.0, 10.0, 10.0, 10.0);
plot(output(:,1), output(:,2)/1e6, 'r');
figure;
plot(targetData(:,1), targetData(:,2)/1e6);

%plot(cal.Energy(:,1), cal.Energy(:,2), 'k')
%hold on;
%plot(energy(:,1), energy(:,2), 'r')