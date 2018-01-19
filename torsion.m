% Torsion Analysis

clear all
close all
clc

Mr = load('torsion_redesign.csv');
Ma = load('torsion_aspire.csv');

spanr = Mr(:,1);
spana = Ma(:,1);
Cmr = Mr(:,2);
Cma = Ma(:,2);
q = 10.04;

cr = 9.4-5.64/49.35*abs(spanr);
ca = zeros(40,1);
for i=1:40
    if abs(spana(i)) <= 20
        ca(i) = 8;
    elseif abs(spana(i)) >= 8
        ca(i) = 8-1.9/18*(abs(spana(i))-20);
    end
end

Tr = q*Cmr.*cr.^2;
Ta = q*Cma.*ca.^2;

plot(spanr,Tr,'b')
hold on
plot(spana,Ta,'r')