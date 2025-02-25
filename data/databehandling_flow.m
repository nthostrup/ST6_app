%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: C:\Users\Bruger\OneDrive\Dokumenter\Uni\Sundhedsteknologi\6. semester\Projekt\data\data_for_mathias_v1.txt
%
% Auto-generated by MATLAB on 21-Apr-2020 15:53:45

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["FlowLs", "Var2", "Var3", "Var4"];
opts.SelectedVariableNames = "FlowLs";
opts.VariableTypes = ["double", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var2", "Var3", "Var4"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var2", "Var3", "Var4"], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable("data/data_for_mathias_v1.txt", opts);

% Convert to output type
FlowLs = tbl.FlowLs;

% Clear temporary variables
clear opts tbl
%% 
%Data from DAN:

y = FlowLs;
x = linspace(0,size(y,1)/2000,size(y,1));
plot(x,y)
xlim([0 18]);
title("Plot af r� data");
xlabel("Tid(s)")
ylabel("Flow (L/s)")
legend("Flow, r� data")
grid on;

%% filtrering
sF = 2000; %sample frequency
Fc=12; % Cut-off frequency (from 2 Hz to 6 Hz depending to the type of your electrod)
N=4; % Filter Order
Wn=2*Fc/sF;
[B, A] = butter(N,Wn, 'low'); %filter's parameters 

yFilt=filtfilt(B,A,y); 
figure;
hold on;
plot(x,yFilt)
xlim([0 18]);
title("Plot af filtreret data (12Hz, 4. orden, lavpas)");
xlabel("Tid(s)")
ylabel("Flow (L/s)")
legend("Flow, filtreret data")
grid on;

%% Calculations
figure;
findpeaks(yFilt,'MinPeakDistance', 3*sF, 'MinPeakHeight', 0.25);
[pks, loc]=findpeaks(yFilt,'MinPeakDistance', 3*sF, 'MinPeakHeight', 0.25);

RR = 60/(x(end)/length(pks));
avgRespDuration = x(end)/length(pks);

%Initializing algorithm values
numberOfBreathsToFind = 125;
endIn = 1;
ends = zeros(1,numberOfBreathsToFind);
starts = zeros(1,numberOfBreathsToFind);
dursIn = zeros(1,numberOfBreathsToFind);
dursExp = zeros(1,numberOfBreathsToFind);
VolsIn = zeros(1,numberOfBreathsToFind);
VolsExp = zeros(1,numberOfBreathsToFind);

for i = 1:numberOfBreathsToFind % Custom amount, based on number of breaths recorded + some to account for errors..
    if endIn == 1 % To account for the first scenario where endIn=1 to NOT get offset of 1
        startIn = find(yFilt(endIn:end) > 0, 1);
        starts(i)=startIn;
    else
        startIn = find(yFilt(endIn:end) > 0, 1) + endIn;
        if(startIn > (endIn+1.5*sF))%Threshold between end of breath to finding start.
            starts(i)=startIn;
            %Finding duration of expiration
            durEx = (startIn-endIn)/sF;
            dursExp(i)=durEx;
            %Finding volume of expiration
            V_exp = trapz(yFilt(endIn:startIn))/sF;
            VolsExp(i)=V_exp;
%         else % Tried to account for error in detection, when 
%            if(i>1 && ~isempty(startIn))
%             ends(i-1) = startIn;
%            end
        end
        
    end
endIn = find(yFilt(startIn:end) < 0, 1) + startIn; % Can be improved by setting threshold between start and end with if.
if(endIn > (startIn+1.5*sF))%Threshold between start of breath to finding end.
    ends(i)=endIn;
    %Finding duration of inspiration
    dur = (endIn-startIn)/sF; % TODO replace StartIn with start(i-1),
    dursIn(i)=dur;

    %Finding volume of inspiration
    V_in = trapz(yFilt(startIn:endIn))/sF;
    VolsIn(i)=V_in;
end
end

starts = starts(find(starts > 0));
ends= ends(find(ends > 0));
dursIn = dursIn(find(dursIn > 0));
dursExp = dursExp(find(dursExp > 0));
VolsIn = VolsIn(find(VolsIn > 0));
VolsExp = VolsExp(find(VolsExp < 0));

figure;
plot(yFilt)
hold on;
plot(starts,yFilt(starts),'*')
plot(ends,yFilt(ends),'*')
legend("signal","starts","ends")

%% Calulations for report
meanDursIn = mean(dursIn);
stdDursIn = std(dursIn);
minDursIn = min(dursIn);
maxDursIn = max(dursIn);

stdDursExp = std(dursExp);
meanDursExp = mean(dursExp);
minDursExp = min(dursExp);
maxDursExp = max(dursExp);

meanVolsIn = mean(VolsIn);
stdVolsIn = std(VolsIn);
minVolsIn = min(VolsIn);
maxVolsIn = max(VolsIn);

meanVolsExp = mean(VolsExp);
stdVolsExp = std(VolsExp);
minVolsExp = min(VolsExp);
maxVolsExp = max(VolsExp);

IE_rel = dursIn(1:end-1)./dursExp;

%Visualize all inspirations:
% for i=1:length(starts)
%     difference = ends(i)-starts(i)+1;
%     x = linspace(0,difference/2000,difference);
%     
%     plot(x,yFilt(starts(i):ends(i)))
%     hold on;
% end

%calculation of mean flow
sum=zeros(1,length(starts));
for i=1:length(starts)
    sum(i)= mean(yFilt(starts(i):ends(i)));
end
Qv_avg = mean(sum);

%% Bodeplot/filterplot
[H, W] = freqz(B,A,1000);
figure;
plotHandle = semilogx((sF*W)/(2*pi),abs(H));
xlim([0 70]);
title("Bodeplot af filter (12Hz, 4. orden, lavpas)");
xlabel("Frekvens (Hz)")
ylabel("Gain")
grid on;
%text(1,1.05,"1Hz")
%line([1 1],[0 1.2],'color','k','Linestyle','--')
%% frekvensplot af r� signal
figure;
fy = abs(fft(detrend(y)));
xf = linspace(0,sF,length(fy));
subplot(1,2,1)
stem(xf,fy)
title("Frekvensplot af r� data");
xlabel("Frekvens (Hz)")
ylabel("Power")
legend("Frekvensplot")
grid on;

subplot(1,2,2)
stem(xf,fy)
xlim([0 70]);
title("Frekvensplot af r� data [0-70Hz]");
xlabel("Frekvens (Hz)")
ylabel("Power")
legend("Frekvensplot")
grid on;
