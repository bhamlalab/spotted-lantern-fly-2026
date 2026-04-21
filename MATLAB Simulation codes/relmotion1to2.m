function [output,time] = relmotion1to2(input1,w1,input2,w2,numcycles,numpts)
wmin=max(w1,w2);
output=[];
f=wmin;%rad/s to 1/s
T=1/f;%Get the period
tt=linspace(0,numcycles*T,numpts);
time=tt;
for i=1:length(tt)
    output=[output input1(tt(i))-input2(tt(i))];
end

end