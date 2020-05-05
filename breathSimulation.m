function [Qv,P] = breathSimulation(RR,Compl, MIP, breathDuration, IMTmodstand, sampleFrekvens)
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
V0_exp = pMax/E;
V0_ins = 0;


t = linspace(0,Ttot,nSamples);
x = t;
Pmus = zeros(1,length(t));
Qv = zeros(1,length(t));
V =  zeros(1,length(t));

dt=1/sF;


for i = 1:length(t)
            if t(i) <= Ti_E*Ttot
                Pmus(i) = pMax*(1-exp(-(1/tauC)*t(i)));
                Qv(i) = Pmus(i)/R_in;
                V(i) = Pmus(i)/E;
                
                ti_end = t(i);
                index_i_end = i;
                Pmus_end = Pmus(i);
                
                V0_exp = V(i);
            
            elseif Ti_E*Ttot < t(i) && t(i) <= Ttot
                Pmus(i) = Pmus_end*(exp(-(1/tauR)*(t(i)-ti_end)));  %(t(i)-ti_end) to make expiration start as if from t=0
                Qv(i) = -1*Pmus(i)/R_in; %Expiratory flow, by definition negative
                V(i) = V(i-1)+Qv(i)*dt;
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




