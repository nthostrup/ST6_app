%% T1.1 randomiseret test af den fysiologiske model
amountToGenerate = 10;
R = randi([16 33],[1 amountToGenerate]);
C = randi([43 75],[1 amountToGenerate])./1000;
pmax = randi([49 105],[1 amountToGenerate]);

for i=1:amountToGenerate
    [Q Pmus] = breathSimulation(R(i),C(i),pmax(i),10,0,1000);
    for j = 1:10 %times to compare
        pmusCalc(i,j)=pmus2(j,R(i),C(i),pmax(i));
        pmusGen(i,j) = Pmus(j*1000);
    end
end

plot(pmusCalc(2,:)); hold on; plot(pmusGen(2,:))

%% T.2.1: Test af detektering af åndedræt case 1 ægte signal fra Dan
%Testinputopsætning
sF = 2000;
data = importdata("data/data_for_mathias_v1.txt",'\t',1);
data = data.data(1:60*sF,1);
    Fc=12; % Cut-off frequency (from 2 Hz to 6 Hz depending to the type of your electrod)
    N=4; % Filter Order
    Wn=2*Fc/sF;
    [B, A] = butter(N,Wn, 'low'); %filter's parameters 

    data=filtfilt(B,A,data); 
startIn = 0;
endIn = 1;
starts = zeros(1,10);
ends = zeros(1,10);
    
detectionCounter = 1;

%Forventetoutput
EXP_OUT_startIn = 112928;
EXP_OUT_endIn = 116673;
EXP_OUT_starts = [878, 12791, 28187, 43072, 54247, 67688, 80036, 99215, 112176];
EXP_OUT_ends=[4584, 18875, 32722, 46190, 58717, 72577, 84645, 103502, 116673]; 
%EXP_OUT_data; 
EXP_OUT_detectionCounter = 9;
j = 1;
k = 1;
%Testen
stepSize = 100;
for i=1:stepSize:length(data)
    [tempStartIn, tempEndIn, tempStarts, tempEnds, detectionCounter] = UNIT_test_parameterDetection(startIn, endIn, starts, ends, data(1:i), detectionCounter, sF);
    if(~(tempStartIn == startIn))
        startIn = tempStartIn;
        starts(j)=startIn;
        j=j+1;
    end
    if(~(tempEndIn == endIn))
        endIn = tempEndIn;
        ends(k)=endIn;
        k=k+1;
    end
    
    clear tempStarts tempEnds;
end
detectionCounter = detectionCounter -1;


%% T.2.1: Test af detektering af åndedræt case 2 signal fra model
%dataopsætning
sF = 1000;
amount = 10;
duration = 10;

R = randi([16 33],[1 amount]);
C = randi([43 75],[1 amount])./1000;
pmax = randi([49 105],[1 amount]);
%% 
close all;
%data setup
flows = zeros(amount,sF*duration); %Each breath with flow simulation gets its own column in the matrix 
pressures = zeros(amount,sF*duration);
numberToRun = 1;
%Generate breaths
for i = 1:amount
    [flows(i,:), pressures(i,:)] = breathSimulation(R(numberToRun), C(numberToRun), pmax(numberToRun), duration,0, sF);
end

%Concatenate breaths to one long sample
data = [];
for i=1:amount
    data = [data flows(i,:)];
end
%Testinputopsætning
startIn = 0;
endIn = 1;
starts = zeros(1,amount);
ends = zeros(1,amount);
    
detectionCounter = 1;

%Forventetoutput
EXP_OUT_starts = [1 10002 20007:10000:100000];
%EXP_OUT_ends= duration*sF*(1/4):duration*sF:duration*amount*sF; 
EXP_OUT_ends=[2940 12945:10000:100000] 
EXP_OUT_detectionCounter = amount;
j = 1;
k = 1;
%Testen
stepSize = 100;
for i=1:stepSize:length(data)
    [tempStartIn, tempEndIn, tempStarts, tempEnds, detectionCounter] = UNIT_test_parameterDetection(startIn, endIn, starts, ends, data(1:i), detectionCounter, sF);
    if(~(tempStartIn == startIn))
        startIn = tempStartIn;
        starts(j)=startIn;
        j=j+1;
    end
    if(~(tempEndIn == endIn))
        endIn = tempEndIn;
        ends(k)=endIn;
        k=k+1;
    end
    
    clear tempStarts tempEnds;
end
detectionCounter = detectionCounter -1;

%Statistik
[Hs Ps] = ttest2(starts, EXP_OUT_starts);
[He Pe] = ttest2(ends,EXP_OUT_ends);
EXP_OUT_starts'
starts'
EXP_OUT_ends'
ends'
plot(data)
grid on;