close all
clear all
clc

codeFolder = '/Users/akira/Documents/Github/Afferented-Muscles-Testing';
dataFolder = '/Volumes/DATA2/TwoMuscleSystemData/CombinedModel/ModelTesting';

modelParameter.pennationAngle = 10*pi/180; %9.6
modelParameter.optimalLength = 3; % 6.8
modelParameter.tendonSlackLength = 0.5; % 24.1
modelParameter.mass = 0.001;
modelParameter.muscleInitialLength = 3; % muscle initial length
modelParameter.tendonInitialLength = 0.5;
load('alpha')
modelParameter.alpha = alpha;
simulationParameter.Lce = 1;

Fs = 10000;
t = 0:1/Fs:5;
amp_temp = 0.1:0.1:1;
forceLevels = 10:10:100;

for i = 1 %:length(amp_temp)
amp = 0.2; %amp_temp(i);
trialN = i;
input = [zeros(1,1*Fs) amp*[0:1/Fs:1] amp*ones(1,length(t)-1*Fs-length(amp*[0:1/Fs:1]))];

output_Song = muscleModel_Song_v2(t,Fs,input,modelParameter);


% fileName = ['Force_' num2str(forceLevels(i))];
% load (fileName)
% Force_MSMS = simout.signals.values;
% 
meanForce_Song(i) = mean(output_Song.Force_tendon(4*Fs:end));
% meanForce_MSMS(i) = mean(Force_MSMS(50000:end));
%%
figure(1)
plot(t,output_Song.Force_tendon)
title('Tendon Force')

figure(2)
plot(t,output_Song.Force_total)
title('Muscle Force')

% figure(3)
% subplot(1,2,1)
% plot(t,output_Song.Lce)
% subplot(1,2,2)
% plot(simout2.signals.values)
% title('Muscle Length')

end

% figure(5)
% plot(forceLevels,meanForce_Song)
% hold on 
% plot(forceLevels,meanForce_MSMS)
% 


