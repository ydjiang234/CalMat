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
backbone(:,2) = backbone(:,2) * 1.0e6;
targetData = load('TargetData.txt');%Displacement (mm) - Load (kN) 
targetData(:,1) = targetData(:,1) / L;
targetData(:,2) = targetData(:,2) * L * 1.0E3;
cal = CalMat( A, I, L, Naxial, revK, backbone, targetData, ampFactor);

filePath = 'ModifiedCPTemplate.tcl';
outPath = 'rendered.tcl';
str = cal.renderTemplate(filePath, outPath, 1.0, 1.0, 1.0, 1.0);
%cal.runOpenSees(filePath)
%out = cal.convertDisp();