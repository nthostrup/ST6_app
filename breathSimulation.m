function [Qv,P] = breathSimulation(RR,Compl, MIP, breathDuration, IMTmodstand, sampleFrekvens)
%breathSimulation: Function to model pressure and flow of breath
%
modst_threshold = IMTmodstand;
R_in = RR;
C = Compl;
pMax = MIP;
tauC = R_in*C;
tauR = tauC;
Ttot = breathDuration;
Ti_E = 1/4;
Ti = Ttot*Ti_E;
Te = Ttot-Ti;
sF = sampleFrekvens;%Static
nSamples = Ttot * sF;


t = linspace(0,Ttot,nSamples);
x = t;
pmusRes = zeros(1,length(t));

for i = 1:length(t)
            if t(i) <= Ti_E*Ttot
                Pmus(i) = pMax*(1-exp(-(1/tauC)*t(i)));
                
                ti_end = t(i);
                index_i_end = i;
                Pmus_end = Pmus(i);
            
            elseif Ti_E*Ttot < t(i) && t(i) <= Ttot
                Pmus(i) = Pmus_end*(exp(-(1/tauR)*(t(i)-ti_end)));  %(t(i)-ti_end) to make expiration start as if from t=0
            end
end

Qv = (Pmus/R_in);
P = Pmus;
end

