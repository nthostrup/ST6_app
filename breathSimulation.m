function [Qv,P, V] = breathSimulation(RR,Compl, MIP, breathDuration, IMTmodstand, sampleFrekvens)
%breathSimulation: Function to model pressure and flow of breath
%
modst_threshold = IMTmodstand; %TODO impl. simulering med IMT threshold
R_in = RR;
C = Compl;
E = 1/C;
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
Pmus = zeros(1,length(t));
Qv = zeros(1,length(t));
V =  zeros(1,length(t));

dt=1/sF;

%Calculating initial values
Pmus(1) = pMax*(1-exp(-(1/tauC)*t(1)));
V(1) = 0; 
Qv(1) = (Pmus(1)-E*V(1))/R_in;

for i = 2:length(t)
            if t(i) <= Ti_E*Ttot
                Pmus(i) = pMax*(1-exp(-(1/tauC)*t(i)));
                V(i) = V(i-1) + Qv(i-1)*dt; 
                Qv(i) = (Pmus(i)-E*V(i))/R_in;
                                
                ti_end = t(i);
                Pmus_end = Pmus(i);
                            
            elseif Ti_E*Ttot < t(i) && t(i) <= Ttot
                Pmus(i) = Pmus_end*(exp(-(1/tauR)*(t(i)-ti_end)));  %(t(i)-ti_end) to make expiration start as if from t=0
                V(i) = V(i-1) + Qv(i-1)*dt; 
                Qv(i) = (Pmus(i)-E*V(i))/R_in;
            end
end

dt = 1/sF;                   % seconds per sample
t2 = (0:dt:Ttot-dt)';     % seconds
%%Sine wave:
Fc = 50;                     % hertz
A = 0.1; % amplitude of noise
noise = A*sin(2*pi*Fc*t2);


Qv = Qv + noise';
P = Pmus;
end




