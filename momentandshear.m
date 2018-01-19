% Bending Moment and Shear Force
% Glider Redesign
% 4 Degrees Angle of Attack

clear all
close all
clc

BMr = load('bendingmoment_redesign.csv');
BMa = load('bendingmoment_aspire.csv');
spanr = BMr(:,1);
Mr = BMr(:,2);
spana = BMa(:,1);
Ma = BMa(:,2);

hold on
plot(spanr,Mr,'b');
plot(spana,Ma,'r');

fMr = spline(spanr,Mr);
fVr = fnder(fMr);
fMa = spline(spana,Ma);
fVa = fnder(fMa);

Vr = abs(ppval(spanr,fVr));
Va = abs(ppval(spana,fVa));

figure
plot(spanr,Vr,'b')
hold on
plot(spana,Va,'r')