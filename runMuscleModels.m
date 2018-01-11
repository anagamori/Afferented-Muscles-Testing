% close all
% clear all
% clc

codeFolder = '/Users/akira/Documents/Github/Afferented-Muscles-Testing';
dataFolder = '/Volumes/DATA2/TwoMuscleSystemData/CombinedModel/ModelTesting';

modelParameter.pennationAngle = 9.6*pi/180;
modelParameter.optimalLength = 6.8;
modelParameter.tendonSlackLength = 24.1;
modelParameter.mass = 0.15;
modelParameter.muscleInitialLength = 6.8; % muscle initial length
modelParameter.tendonInitialLength = 24.1;
load('alpha')
modelParameter.alpha = alpha;
simulationParameter.Lce = 1;

Fs = 20000;
t = 0:1/Fs:5;
amp_temp = 0.1:0.1:1;

for i = 10 %length(amp_temp)
amp = amp_temp(i);
trialN = i;
input = [zeros(1,1*Fs) amp*[0:1/Fs:1] amp*ones(1,3*Fs)];
%input = 0.3*sin(2*pi*1.*t-pi/2)+0.4;
    
output_1 = muscleModel_Song(t,Fs,input,modelParameter);
% output_2 = muscleModel_Tsianos(t,Fs,input,modelParameter,simulationParameter);
output_3 = muscleModel_Combined(t,Fs,input,modelParameter,1);

maxForce1(i) =  mean(output_1.Force_tendon(4*Fs:5*Fs));
% maxForce2(i) = mean(output_2.Force_total(6*Fs:7*Fs));
maxForce3(i) = mean(output_3.Force_total(4*Fs:5*Fs));
CoV(i) = std(output_3.Force_total(4*Fs:5*Fs))/mean(output_3.Force_total(4*Fs:5*Fs));

%%
figure(1)
plot(t,output_1.Force_tendon)
hold on 
% plot(t,output_2.Force_total)
% hold on 
plot(t,output_3.Force_total)
hold on 
%plot(t,output_3.Force_total)

cd (dataFolder)
save(['output_Song_' num2str(trialN)],'output_1')
%save(['output_Tsianos_' num2str(trialN)],'output_2')
% save(['output_MU_' num2str(trialN)],'output_3','-v7.3')
cd (codeFolder)
end

%figure(2)
% plot(amp_temp,maxForce1./maxForce1(end))
% hold on
% plot(amp_temp,maxForce2./maxForce2(end))
% hold on
%plot(amp_temp,maxForce3./maxForce3(end))
%legend('Song','Tsianos','Motor Unit')

% figure(4)
% plot(amp_temp,CoV)

