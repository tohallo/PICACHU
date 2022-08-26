%% PICArroCHamberflUx (PICACHU) v0.2
%PICACHU.m is the main program.  In the current configuration, this file,
%and the other .m files (functions) it refers to, should be copied
%alongside the data folders it refers to, and under a master sample-date
%based folder.  this way, the PICACHU.m file can be customized for each
%date as needed.  in the future, it may be possible to streamline this so
%only data files need to change

%Based on the LICOR flux equations from here
%https://www.licor.com/env/support/LI-8100A/topics/deriving-the-flux-equation
%https://www.licor.com/env/support/LI-6800/topics/soil-chamber-theory.html
% /\︿╱\
% \0_ 0 /╱\╱
% \▁︹_/

clear all
close all
%% constants
R = 8.314; % J K-1 mol ‑1
irc = 9.84; %inner radius chamber (cm) per Georgia Seyfried, 8/24/22
S = pi*(irc^2); %exposed soil area (cm2)
offsets = [7]; %chamber base offsets (cm)
hc = 91.44; %nominal chamber height (cm) per Georgia Seyfried, 8/24/22 (3 feet)
v8100 = 19;  %(cm^3)
vPicarro = 55; %(cm^3)
vtubing = 14.7; %for 100' of 1/8" ID tubing
V = S.*(offsets+hc) + vPicarro + v8100 + vtubing; %total chamber volume including offset (cm^3)  (volume of a cylinder plus known system volumes)
P = 101.325; % is atmospheric pressure (kPa) (get from tower)
W = 30; %water mole fraction in the chamber (mmol mol-1) (get from chamber RH)
Tcham = 35; %chamber air temperature (°C) (get from chamber)
% C' is dry CO2 concentration (mol mol -1),

%% read Picarro dats
%d = readtable('.\2022-06-14-HB4\JFAADS2140-20220614-093601-DataLog_User.dat');
tt= combine_Picarro_DAT(pwd,'Picarro',datenum('2022-06-13'),datenum('2022-06-15')); %a function Tom wrote
%tt = table2timetable(d, 'RowTimes', datetime((d.DATE) + duration(d.TIME),'format','yyyy-MM-dd HH:mm:ss.SSS'));

% figure(1)
% plot(tt.Time,tt.CH4_dry,'k.-')

%% read Campbell (CSI) logger

csi = readtimetable('Campbell-logger/CR800_Table1_5sec.dat');
csi.TIMESTAMP = csi.TIMESTAMP+hour(1); %need to adjust clocks??

%% read metadata (time stamps and chamber IDs)

meta = readtable('metadata/Picarro-metadata-2022-06-14.csv');
m = table2timetable(meta, 'RowTimes', datetime(meta.date,'convertfrom','yyyyMMdd') + duration(meta.time_start));% convert to a time table where Time will be the chamber start time
m.endtime = datetime(meta.date,'convertfrom','yyyyMMdd') + duration(meta.time_end); %create a variable in Matlab datetime format for the chamber end time

isplant = contains(m.chamberID,'plant'); %logical index true where 'plant' occurs in chamber ID

%loop through each line of metadata
for i = 1:height(m)
    inds = timerange(m.Time(i)+minutes(1.25),m.endtime(i)-minutes(0.5)); %an index of rows that are between the start (inclusive) and end (exclusive) times
    %     figure(1)
    %     hold on
    %     plot(tt.Time(inds),tt.CH4_dry(inds),'ro') %update the plot with red dots for measurement times
    %     hold off
    dumtime = datenum(datevec(tt.Time(inds))); %convert meas. time to matlab serial date format (fractional days)
    x = dumtime-floor(dumtime); %subtract off the integer value of the date to get fractional day
  
    %% CO2 FLUX section
    % y = tt.CH4_dry(inds); %extract CH4 meas.
    y = tt.CO2_dry(inds); %extract CO2 meas.
    
    figure(i+1)
    set(gcf,'position',[1100 870 1100 468])
    subplot(1,2,1)
    fid = gcf;
    scatter(x,y,'ro')
    stats = placestatsonfigure_robust(x,y,1,'k',fid); % a function Tom wrote to fit linear regression and plot it.  robust ignores outliers
    title(m.chamberID(i))
    xlabel('fractional day')
    ylabel(['CO_2'])
    
    %%fit exponential model
    xsecs = (x-x(1))*24*60*60; %elapsed seconds
    
    inds = timerange(m.Time(i)+minutes(1.25),m.endtime(i)-minutes(0.5)); %an index of rows that are between the start (inclusive) and end (exclusive) times
    
    %   Beta = nlinfit(X,Y,@(beta,X)modelfun(beta, c, X),beta0)
    X = [xsecs];
    c0 = y(1); %initial CO2, fit linear in future
    %beta0 = [0.002 2.9]; %methane
    if stats.slope < 0
        beta0 = [0.003 350]; %CO2, negative slope
    else
        beta0 = [1e-6 600];
    end
    
    options = statset('Robust','off','MaxIter',100000,'Display','final');
    [beta1,r1,J1,COVB1,mse1] = nlinfit(xsecs,y,@(beta,X)chamber_exp2(beta,X,c0,0),beta0,options);
    %[C,delta] = nlpredci(@chamber_exp2,xsecs,beta1,r1,'covar',COVB1);
    
    C = chamber_exp2(beta1,xsecs,y(1),0);
    hold on
    plot(x,C,'g-')
    
    params.co2(i,:) = beta1;
    
    a = params.co2(i,1);
    cx = params.co2(i,2);
    dcdt  = a.*(cx-c0); %eq 11-15
    
    if isplant(i)
        vol = V; %plant chamber volume (needs offsets!)
    else
        vol = 10315; %floating chamber volume (cu cm); per Georgia Seyfried on 8/24/22
    end
    
    F = (10.*vol.*P.*(1-(W/1000)).* dcdt) ./ (R.*S.*(Tcham+273.15)); %main flux equation
    fluxes.co2(i,1) = F; %CO2 in column 1
    
    %% CH4 section
      y = tt.CH4_dry(inds); %extract CH4 meas.
    %y = tt.CO2_dry(inds); %extract CH4 meas.
    
    %figure(i+1)
    subplot(1,2,2)
    fid = gcf;
    scatter(x,y,'bo')
    stats = placestatsonfigure_robust(x,y,1,'k',fid); % a function Tom wrote to fit linear regression and plot it.  robust ignores outliers
    title(m.chamberID(i))
    xlabel('fractional day')
    ylabel(['CH_4'])
    
    %%fit exponential model
    xsecs = (x-x(1))*24*60*60; %elapsed seconds
    
    inds = timerange(m.Time(i)+minutes(1.25),m.endtime(i)-minutes(0.5)); %an index of rows that are between the start (inclusive) and end (exclusive) times
    
    %   Beta = nlinfit(X,Y,@(beta,X)modelfun(beta, c, X),beta0)
    X = [xsecs];
    c0 = y(1); %initial CO2, fit linear in future
    beta0 = [0.002 2.9]; %methane
%     if stats.slope < 0
%         beta0 = [0.003 350]; %CO2, negative slope
%     else
%         beta0 = [1e-6 600];
%     end
    
    options = statset('Robust','off','MaxIter',100000,'Display','final');
    [beta1,r1,J1,COVB1,mse1] = nlinfit(xsecs,y,@(beta,X)chamber_exp2(beta,X,c0,0),beta0,options);
    %[C,delta] = nlpredci(@chamber_exp2,xsecs,beta1,r1,'covar',COVB1);
    
    C = chamber_exp2(beta1,xsecs,y(1),0);
    hold on
    plot(x,C,'g-')
    
    params.ch4(i,:) = beta1;
    
    a = params.ch4(i,1);
    cx = params.ch4(i,2);
    dcdt  = a.*(cx-c0); %eq 11-15
    
    if isplant(i)
        vol = V; %plant chamber volume (needs offsets!)
    else
        vol = 10315; %floating chamber volume (cu cm); per Georgia Seyfried on 8/24/22
    end
    
    F = (10.*vol.*P.*(1-(W/1000)).* dcdt) ./ (R.*S.*(Tcham+273.15));
    fluxes.ch4(i,1) = F; 
    
    %% figures
    
    figure(i+48)
    subplot(1,3,1)
    plot(csi.TIMESTAMP(inds),[csi.TairfromRH_Avg(inds) csi.Tair_C_Avg(inds)],'o-')
    ylabel('Tair chamber (oC)')
    subplot(1,3,2)
    plot(csi.TIMESTAMP(inds),csi.RH_Avg(inds),'o-')
    ylabel('RH (%)')
    subplot(1,3,3)
    plot(csi.TIMESTAMP(inds),csi.diffPress_mv_Avg(inds),'o-')
    ylabel('diff pressure')
    
    %export_fig('Annandale-06-14-2022v2.pdf','-append','-painters')
    
end

%% 
figure('position',[1100 870 1100 468])
subplot(1,2,1)
bar(fluxes.co2)
grid on
title('Annandale - June 14, 2022')
xlabel('measurement #')
ylabel('CO_2 flux (\mumol CO_2 m^-^2 sec^-^1)')

subplot(1,2,2)
bar(fluxes.ch4)
grid on
title('Annandale - June 14, 2022')
xlabel('measurement #')
ylabel('CH_4 flux (\mumol CH_4 m^-^2 sec^-^1)')



% ⠸⣷⣦⠤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣠⣤⠀⠀⠀
% ⠀⠙⣿⡄⠈⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠔⠊⠉⣿⡿⠁⠀⠀⠀
% ⠀⠀⠈⠣⡀⠀⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠊⠁⠀⠀⣰⠟⠀⠀⠀⣀⣀
% ⠀⠀⠀⠀⠈⠢⣄⠀⡈⠒⠊⠉⠁⠀⠈⠉⠑⠚⠀⠀⣀⠔⢊⣠⠤⠒⠊⠉⠀⡜
% ⠀⠀⠀⠀⠀⠀⠀⡽⠁⠀          ⠀⠩⡔⠊⠁⠀⠀⠀⠀⠀⠀⠇
% ⠀⠀⠀⠀⠀⠀⠀⡇⢠⡤⢄⠀⠀⠀⠀⠀⡠⢤⣄⠀⡇⠀⠀⠀⠀⠀⠀⠀⢰⠀
% ⠀⠀⠀⠀⠀⠀⢀⠇⠹⠿⠟⠀⠀⠤⠀⠀⠻⠿⠟⠀⣇⠀⠀⡀⠠⠄⠒⠊⠁⠀
% ⠀⠀⠀⠀⠀⠀⢸⣿⣿⡆⠀⠰⠤⠖⠦⠴⠀⢀⣶⣿⣿⠀⠙⢄⠀⠀⠀⠀⠀⠀
% ⠀⠀⠀⠀⠀⠀⠀⢻⣿⠃⠀⠀⠀⠀⠀⠀⠀⠈⠿⡿⠛⢄⠀⠀⠱⣄⠀⠀⠀⠀
% ⠀⠀⠀⠀⠀⠀⠀⢸⠈⠓⠦⠀⣀⣀⣀⠀⡠⠴⠊⠹⡞⣁⠤⠒⠉⠀⠀⠀⠀⠀
% ⠀⠀⠀⠀⠀⠀⣠⠃⠀⠀⠀⠀⡌⠉⠉⡤⠀⠀⠀⠀⢻⠿⠆⠀⠀⠀⠀⠀⠀⠀
% ⠀⠀⠀⠀⠀⠰⠁⡀⠀⠀⠀⠀⢸⠀⢰⠃⠀⠀⠀⢠⠀⢣⠀⠀⠀⠀⠀⠀⠀⠀
% ⠀⠀⠀⢶⣗⠧⡀⢳⠀⠀⠀⠀⢸⣀⣸⠀⠀⠀⢀⡜⠀⣸⢤⣶⠀⠀⠀⠀⠀⠀
% ⠀⠀⠀⠈⠻⣿⣦⣈⣧⡀⠀⠀⢸⣿⣿⠀⠀⢀⣼⡀⣨⣿⡿⠁⠀⠀⠀⠀⠀⠀
% ⠀⠀⠀⠀⠀⠈⠻⠿⠿⠓⠄⠤⠘⠉⠙⠤⢀⠾⠿⣿⠟⠋




