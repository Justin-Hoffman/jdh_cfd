clear all;clc;close all;
t=-pi:.1:pi;
for i=1:length(t)
    Alison(i)=16*(sin(t(i)))^3;
    Justin(i)=13*cos(t(i))-5*cos(2*t(i))-2*cos(3*t(i))-cos(4*t(i));
end
plot(Alison,Justin,'r'); xlabel('Alison'); ylabel('Justin'); title('Meant To Be');