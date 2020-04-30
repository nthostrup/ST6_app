%% fV = 10, P01 = 0.5: Fresnel %p01 er trykket målt 100ms efter igangsættelse af inspiration.
close all;
clear
fV = 10;
p01 = 4;

Ttot = 60/fV;
pMax = p01/(1-exp(-(0.1*(fV + 4*p01))/10));
Ti_tot = 0.0125*fV + 0.125;
Ti = Ttot*Ti_tot;
Te = Ttot-Ti;
tauC = 10/(fV + 4*p01);
tauR = 10/(fV + p01/2);

syms t y;
pmus1(t) = -pMax*(1-exp(-(1/tauC)*t)); %0 < t <= Ti
pmus2(t) = exp(-(1/tauR)*t);  %Ti < t <= Ttot
pmus3(t,y) = pmus1(t) + pmus2(y);

ttime = linspace(0,Ttot,200);
for i = 1:length(ttime)
    if(ttime(i) <= Ti)
        pmusRes(i) = pmus1(ttime(i));
        ti_end = ttime(i);
    elseif(ttime(i) > Ti && ttime(i) <= Ttot)
        pmusRes(i) = pmus1(ti_end)*pmus2(ttime(i)-ti_end);
    end
end

t1 = linspace(0,Ti,100);
t2 = linspace(Ti,Ttot,100);

ty1 = [linspace(0,Ti,100) zeros(1,100)];
ty2 = [zeros(1,100) linspace(Ti,Ttot,100)];

% plot(t1,pmus1(t1),'b')
% plot(t2,pmus2(t2),'b')

%plot([ty1(1:100) ty2(101:200)],pmus3(ty1,ty2))
x = linspace(0,Ttot,200);
plot(pmusRes)
hold on;

%% fV = 10, P01 = 0.5: Fresnel %p01 er trykket målt 100ms efter igangsættelse af inspiration.
fV = 10;
p01 = 10;

Ttot = 60/fV;
pMax = p01/(1-exp(-(0.1*(fV + 4*p01))/10));
Ti_tot = 0.0125*fV + 0.125;
Ti = Ttot*Ti_tot;
Te = Ttot-Ti;
tauC = 10/(fV + 4*p01);
tauR = 10/(fV + p01/2);

syms t y;
pmus1(t) = -pMax*(1-exp(-(1/tauC)*t)); %0 < t <= Ti
pmus2(t) = pmus1(Ti)*exp(-(1/tauR)*t);  %Ti < t <= Ttot
pmus3(t,y) = pmus1(t) + pmus2(y);


ttime = linspace(0,Ttot,200);
pmusRes = zeros(1,length(ttime));
for i = 1:length(ttime)
    if(ttime(i) <= Ti)
        pmusRes(i) = -pMax*(1-exp(-(1/tauC)*ttime(i)));
        ti_end = ttime(i);
        p_end = pmusRes(i);
    elseif(ttime(i) > Ti && ttime(i) <= Ttot)
        pmusRes(i) = p_end*exp(-(1/tauR*(ttime(i)-ti_end)));
    end
end

t1 = linspace(0,Ti,100);
t2 = linspace(Ti,Ttot,100);

ty1 = [linspace(0,Ti,100) zeros(1,100)];
ty2 = [zeros(1,100) linspace(Ti,Ttot,100)];

% plot(t1,pmus1(t1),'b')
% plot(t2,pmus2(t2),'b')

%plot([ty1(1:100) ty2(101:200)],pmus3(ty1,ty2))
x = linspace(0,Ttot,200);
plot(x,pmusRes)
hold on;


%% Dans impl.
t=ttime;
for i = 1:length(t)
            if t(i) <= Ti_tot*Ttot
                Pmus(i) = pMax*(1-exp(-(1/tauC)*t(i)));
                ti_end = t(i);
                Pmus_end = Pmus(i);
            elseif Ti_tot*Ttot < t(i) && t(i) <= Ttot
                Pmus(i) = Pmus_end*(exp(-(1/tauR)*(t(i)-ti_end)));  %(t(i)-ti_end) to make expiration start as if from t=0
            end
end
figure;
plot(Pmus)

%% "Den rigtige fysiologiske model"
R_in = 2.5;
R_ex = R_in;
C = 0.2;
pMax = 50;
tauC = R_in*C;
tauR = R_ex*C;
Ttot = 6;
Ti_E = 1/4;
Ti = Ttot*Ti_E;
Te = Ttot-Ti;
sF = 2000;
nSamples = Ttot*sF;


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
%figure;
%plot(x,Pmus)
hold on;

dt = 1/sF;                   % seconds per sample
t2 = (0:dt:Ttot-dt)';     % seconds
%%Sine wave:
Fc = 50;                     % hertz
A = 0.001; %amplitude
noise = A*sin(2*pi*Fc*t2);


Qv = (diff(Pmus)/((R_in+R_ex)/2)); %% Todo, fix size
Qv = [Qv 0] + noise';
%Qv = (Pmus/((R_in+R_ex)/2))./nSamples; % Divide by samples. %%TODO overvej om der skal regnes med Delta P i stedet? altså diff?
V = trapz(Qv(1:index_i_end));
plot(x,Qv)
%Qv_avg = mean(Qv(1:index_i_end));
%V = Qv_avg*Ti;
