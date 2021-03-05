%****************************************************************
%   MiniProjectOne.m
%
%   PROGRAM DESCRIPTION
%   This program calculates and plots net acceleration and HIC
%   (head injury criterion) values during NHTSA frontal impact
%   testing of a Toyota Corolla, Honda CR-V, and a Ford F-150.
%
%   INPUT: None user-provided
%   OUTPUT: Subplot popups and calculated values in .txt file
%
%   WRITTEN BY: Evan Schimmel
%               01/11/2021 14:27
%
%****************************************************************

clc
clear variables
close all

file_number=fopen('MiniProjectOne.txt','w');

%% DATA IMPORT %%

toyota_corolla=xlsread('toyota_corolla_data'); %matrix created from Toyota Corolla data
honda_cr_v=xlsread('honda_cr_v_data'); %matrix created from Honda CR-V data
ford_f_150=xlsread('ford_f_150_data'); %matrix created from Ford F-150 data
[n_rows,n_cols]=size(toyota_corolla); %values for number of rows and columns of Toyota Corolla Data (because all data sets are the same size, values will be used for all three sets)

%% NET ACCELERATION CALCULATIONS %%

for i=1:n_rows %loop for net accelerations and time conversions
    timeplot(i)=toyota_corolla(i,1)*1000; %creating time vector for plotting (convert s to ms)
    netaccel_toyota_corolla(i)=sqrt((toyota_corolla(i,2)^2)+(toyota_corolla(i,3)^2)+(toyota_corolla(i,4)^2)); %calculating net accelerations for Toyota Corolla (m/s^2)
    netaccel_honda_cr_v(i)=sqrt((honda_cr_v(i,2)^2)+(honda_cr_v(i,3)^2)+(honda_cr_v(i,4)^2)); %calculating net accelerations for Honda CR-V (m/s^2)
    netaccel_ford_f_150(i)=sqrt((ford_f_150(i,2)^2)+(ford_f_150(i,3)^2)+(ford_f_150(i,4)^2)); %calculating net accelerations for Ford F-150 (m/s^2)
end

%% SEVERITY VALUE CALCULATIONS %%

offset=0.015/0.0001; %calculating number of data points to be utilized in each time window (total time/time of each data point)
for win=1:n_rows-offset %loop for severity windows
    sum_toyota_corolla=0; %initialize Toyota Corolla sum to zero
    sum_honda_cr_v=0; %initialize Honda CR-V sum to zero
    sum_ford_f_150=0; %initialize Ford F-150 sum to zero
    t1=toyota_corolla(win,1); %calculate the initial time in the window (s)
    t2=toyota_corolla(win+offset,1); %calculate the final time in the window (s)
    for index=win:win+(offset-1) %loop for trapezoidal integration of accelerations in m/s^2
        sum_toyota_corolla=sum_toyota_corolla+0.0001*((netaccel_toyota_corolla(index)+netaccel_toyota_corolla(index+1))/2); %recursively defined sum for Toyota Corolla integral (m/s)
        sum_honda_cr_v=sum_honda_cr_v+0.0001*((netaccel_honda_cr_v(index)+netaccel_honda_cr_v(index+1))/2); %recursively defined sum for Honda CR-V integral (m/s)
        sum_ford_f_150=sum_ford_f_150+0.0001*((netaccel_ford_f_150(index)+netaccel_ford_f_150(index+1))/2); %recursively defined sum for Ford F-150 integral (m/s)
    end
    severity_toyota_corolla(win)=((((1/(t2-t1))*sum_toyota_corolla)^2.5)*(t2-t1)); %calculate Toyota Corolla severity for given time window
    severity_honda_cr_v(win)=((((1/(t2-t1))*sum_honda_cr_v)^2.5)*(t2-t1)); %calculate Honda CR-V severity for given time window
    severity_ford_f_150(win)=((((1/(t2-t1))*sum_ford_f_150)^2.5)*(t2-t1)); %calculate Ford F-150 severity for given time window
    sevwindow(win)=win; %create a vector to define the numbered time window
end

%% MAXIMUM NET ACCELERATION & HIC VALUE INDEX/VALUE CALCULATIONS %%

netaccelmax_toyota_corolla=max(netaccel_toyota_corolla); %calculate the max net acceleration value of the Toyota Corolla
netaccelmax_honda_cr_v=max(netaccel_honda_cr_v); %calculate the max net acceleration value of the Honda CR-V
netaccelmax_ford_f_150=max(netaccel_ford_f_150); %calculate the max net acceleration value of the Ford F-150

[max_severity_toyota_corolla,ind_severity_toyota_corolla]=max(severity_toyota_corolla); %calculate the index and value of the Toyota Corolla HIC value
[max_severity_honda_cr_v,ind_severity_honda_cr_v]=max(severity_honda_cr_v); %calculate the index and value of the Honda CR-V HIC value
[max_severity_ford_f_150,ind_severity_ford_f_150]=max(severity_ford_f_150); %calculate the index and value of the Ford F-150 HIC value

%% PRINTING DATA %%

%print the net acceleration data for each vehicle to the .txt file
fprintf(file_number,'Comparison of the peak resultant head acceleration experienced:\n\n');
fprintf(file_number,'Toyota Corolla: %5.2fg\n\n',netaccelmax_toyota_corolla);
fprintf(file_number,'Honda CR-V: %5.2fg\n\n',netaccelmax_honda_cr_v);
fprintf(file_number,'Ford F-150: %5.2fg\n\n\n',netaccelmax_ford_f_150);

%print the HIC value data for each vehicle to the .txt file
fprintf(file_number,'Comparison of the calculated head injury criterion (HIC) values:\n\n');
fprintf(file_number,'Toyota Corolla: %4.0f\n\n',max_severity_toyota_corolla);
fprintf(file_number,'Honda CR-V: %4.0f\n\n',max_severity_honda_cr_v);
fprintf(file_number,'Ford F-150: %4.0f',max_severity_ford_f_150);

%% PLOTTING DATA %%

%net acceleration grouping
figure(1)
sgtitle('Resultant Head Acceleration during Frontal Impact Test')

%Toyota Corolla net acceleration subplot
subplot(3,1,1)
plot(timeplot,netaccel_toyota_corolla)
xlabel('Time (ms)')
ylabel('Acceleration (g)')
axis([-50 300 0 60])
legend('Toyota Corolla')

%Honda CR-V net acceleration subplot
subplot(3,1,2)
plot(timeplot,netaccel_honda_cr_v)
xlabel('Time (ms)')
ylabel('Acceleration (g)')
axis([-50 300 0 60])
legend('Honda CR-V')

%Ford F-150 net acceleration subplot
subplot(3,1,3)
plot(timeplot,netaccel_ford_f_150)
xlabel('Time (ms)')
ylabel('Acceleration (g)')
axis([-50 300 0 60])
legend('Ford F-150')

%severity grouping
figure(2)
sgtitle('Variation in Severity during Frontal Impact Test')

%Toyota Corolla severity subplot
subplot(3,1,1)
plot(sevwindow,severity_toyota_corolla,ind_severity_toyota_corolla,max_severity_toyota_corolla,'ok')
xlabel('Time window number')
ylabel('Severity')
axis([0 3500 0 250])
legend('Toyota Corolla','HIC value')

%Honda CR-V severity subplot
subplot(3,1,2)
plot(sevwindow,severity_honda_cr_v,ind_severity_honda_cr_v,max_severity_honda_cr_v,'ok')
xlabel('Time window number')
ylabel('Severity')
axis([0 3500 0 250])
legend('Honda CR-V','HIC value')

%Ford F-150 severity subplot
subplot(3,1,3)
plot(sevwindow,severity_ford_f_150,ind_severity_ford_f_150,max_severity_ford_f_150,'ok')
xlabel('Time window number')
ylabel('Severity')
axis([0 3500 0 250])
legend('Ford F-150','HIC value')

%%
fclose(file_number);