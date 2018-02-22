close all
clear all
clc

forceLevels = 10:10:100;
for i = 1:10
    fileName = ['Force_' num2str(forceLevels(i))];
    load (fileName)
    Force = simout.signals.values;
    
    meanForce(i) = mean(Force(50000:end));
end

figure(1)
plot(forceLevels,meanForce./meanForce(end)*100)
xlabel('Activation (%MVC)')
ylabel('%MVC')


