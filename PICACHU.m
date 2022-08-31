%% PICArroCHamberflUx (PICACHU) v0.2
%PICACHU.m is the main program.  In the current configuration, this file,
%and the other .m files (functions) it refers to, should be copied
%alongside the data folders it refers to, and under a master sample-date
%based folder.  this way, the PICACHU.m file can be customized for each
%date as needed.  in the future, it may be possible to streamline this so
%only data files need to change

%programmed by Tom O'Halloran (tohallo@clemson.edu), summer 2022

%Based on the LICOR flux equations from here
%https://www.licor.com/env/support/LI-8100A/topics/deriving-the-flux-equation
%https://www.licor.com/env/support/LI-6800/topics/soil-chamber-theory.html
% /\︿╱\
% \0_ 0 /╱\╱
% \▁︹_/

%download the following folders with input data in the same diretory as this .m file

%Campbell-logger
%metadata
%Picarro

clear all
close all
%% constants
R=8.31446; %Universal gas constant, J/K/mol
Rd = 287; %gas constant for dry air J kg^-1 K^-1
Cpd = 1004.67; %heat capacity of dry air J kg^-1 K^-1
Eps =0.622; %ratio molecular weight water vapor to dry air
Na=6.02214e23; %Avogadro's number
irc = 9.84; %inner radius chamber (cm) per Georgia Seyfried, 8/24/22
S = pi*(irc^2); %exposed soil area (cm2)
%offsets = [7]; %chamber base offsets (cm) (put into metadata file)
hc = 91.44; %nominal chamber height (cm) per Georgia Seyfried, 8/24/22 (3 feet)
v8100 = 19;  %(cm^3)
vPicarro = 55; %(cm^3)
vtubing = 14.7; %for 100' of 1/8" ID tubing
%V = S.*(offsets+hc) + vPicarro + v8100 + vtubing; %total chamber volume including offset (cm^3)  (volume of a cylinder plus known system volumes)
Vmisc = vPicarro + v8100 + vtubing; %total chamber volume including offset (cm^3)  (volume of a cylinder plus known system volumes)
%P = 101.325; % is atmospheric pressure (kPa) (get from tower)
%W = 30; %water mole fraction in the chamber (mmol mol-1) (get from chamber RH)
%Tcham = 35; %chamber air temperature (°C) (get from chamber)
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

%% read flux tower data
fluxfile = dir('towerflux/*.csv'); %should be just one EddyPro output .csv in this folder
giant=readtable(['towerflux/' fluxfile.name],'DateTimeType','text');
DT1=datetime(join([string(giant.date) string(giant.time)]));
flux = table2timetable(giant,'RowTimes',DT1);
clear giant DT1

press = naninterpe(autoclean(flux.air_pressure))./1000; %[kPa] remove bad values and interpolate gaps (not great)
fco2 = flux.co2_flux;
fco2qc = flux.qc_co2_flux;
fch4 = flux.ch4_flux;
fch4qc = flux.qc_ch4_flux;

%filter on range and using QAQC filter
fco2(fco2<-50 | fco2>50 | fco2qc>1)=NaN;
fch4(fch4<-10 | fch4>10 | fch4qc>1)=NaN;

flux.fco2 = fco2;
flux.fch4 = fch4;
flux.press = press;

%% read tower met data
metfile = dir('towermet/*.csv'); %Compile met .xlxs from Mike into one CSV to cover dates needed
giant=readtable(['towermet/' metfile.name],'DateTimeType','text');
sd = csitime2date(giant.YEAR,giant.JULIEN,giant.Time); %matlab serial date
DT1 = datetime(datevec(sd));
met = table2timetable(giant,'RowTimes',DT1);
clear giant DT1

%take out of table to do some filtering.  this is probably dumb but the
%table notation is miserable
Tair = met.AirT_probe_1;
RHair = met.RH_probe_1;
SWin = met.SW_IN;
%filtering
Tair(Tair<0 | Tair>60)=NaN; %clean bad values
RHair(RHair<10 | RHair>104)=NaN;
RHair(RHair>100)=100;
SWin(SWin>2000 | SWin<-10)=NaN; %(clean bad values)
SWin(SWin<1)=0; %reset very small or slightly negative values to zero (at night)

%put back into table
met.Tair = Tair;
met.RH = RHair;
met.SWin = SWin;

%% read metadata (time stamps and chamber IDs)

meta = readtable('metadata/Picarro-metadata-2022-06-14.csv');
m = table2timetable(meta, 'RowTimes', datetime(meta.date,'convertfrom','yyyyMMdd') + duration(meta.time_start));% convert to a time table where Time will be the chamber start time
m.endtime = datetime(meta.date,'convertfrom','yyyyMMdd') + duration(meta.time_end); %create a variable in Matlab datetime format for the chamber end time
m.offset = m.offset.*2.54; % inches to cm
isplant = contains(m.chamberID,'plant'); %logical index true where 'plant' occurs in chamber ID

%% preprocess chamber met, tower flux, and tower met data in loop
%loop through each line of metadata
for i = 1:height(m)
    
    %this is a new matlab datetime function that creates flexible time
    %indices. the length of the vector will change depending on the
    %frequency of the data.  this is really handy but also a little
    %confusing
    
    %adjust chamber open and close times per dead bands. can improve this.
    starttime = m.Time(i)+minutes(1.25); %from metadata
    endtime = m.endtime(i)-minutes(0.5); %from metadata
    meastime(i) = mean([starttime endtime]); % midpoint of measurement
    
    inds = timerange(starttime,endtime); %an index of rows that are between the start (inclusive) and end (exclusive) times
    %ie, inds can be applied to timetable data at any frequency and it will
    %find the appropriate rows

    
    % CHAMBER DATA
    if ~isempty(csi.Tair_C_Avg(inds))
        Tcham1(i,1) = nanmean(csi.Tair_C_Avg(inds)); %oC thermocouple if needed
        Tcham2(i,1) = nanmean(csi.TairfromRH_Avg(inds));
        RHcham(i,1) = nanmean(csi.RH_Avg(inds)); % percent
        diffpress(i,1) = nanmean(csi.diffPress_mv_Avg(inds)); %inH2O?
        dum = csi.TIMESTAMP(inds);
        timecham(i,1) = dum(1); %first value each chamber time
    else
        Tcham1(i,1) = NaN;
        Tcham2(i,1) = NaN;
        RHcham(i,1) = NaN;
        diffpress(i,1) = NaN;
        timecham(i,1) = NaT;
    end
    
    %FLUX DATA
    co2f(i,1) = interp1(flux.Time,flux.fco2,starttime,'nearest'); %interpolate nearest half-hourly flux to start time of chamber meas.
    ch4f(i,1) = interp1(flux.Time,flux.fch4,starttime,'nearest'); %interpolate nearest half-hourly flux to start time of chamber meas.
    P(i,1) = interp1(flux.Time,flux.press,starttime,'nearest'); %interpolate nearest half-hourly flux to start time of chamber meas.
    h2o(i,1) = interp1(flux.Time,flux.h2o_mole_fraction,starttime,'nearest'); %interpolate nearest half-hourly flux to start time of chamber meas.
    
    % FLUX TOWER MET DATA
    ta(i,1) = nanmean(met.Tair(inds));
    rh(i,1) = nanmean(met.RH(inds));
    swin(i,1) = nanmean(met.SWin(inds));

end

%fill chamber met gaps
Tfinal = Tcham2; %initialize with thermistor/hygristor
Tfinal(~isplant) = Tcham1(~isplant); %replace non-plant (floating) chamber with thermocouple
Tfinal(isnan(Tfinal)) = ta(isnan(Tfinal)); %replace gaps with tower

RHcham(RHcham>100 | RHcham<10) = NaN;
RHfinal = RHcham;
RHfinal(isnan(RHfinal)) = rh(isnan(RHfinal)); %replace missing chamber with tower

% do some calcs with resulting vectors of values

% Saturation vappor pressure ---> es = 611Pa*exp((17.2694*tair)/(tair +237.3))
es = 611.* exp((17.2694.*(Tfinal))./((Tfinal) + 237.3)); % Saturation vaporpressure in Pa (tair in oC)

% Actual vapor pressure ---> e = (RH/100)*es
e = (RHfinal./100).*es; % Actual vapor pressure in Pa

%mixing ratio (k/kg)
w = (Eps*(e./(P - e))).*1000;

figure
plot(meastime,[Tcham1 Tcham2 ta],'o-')
ylabel('chamber air temp (oC)')
legend({'thermocouple','Tair/H','tower'})
title('Chamber properties')

figure
plot(meastime,[RHcham rh RHfinal],'o-')
ylabel('relative humidity  (%)')
legend({'Tair/H','tower','final'})
title('Chamber properties')

figure
plot(meastime,diffpress,'o')
ylabel('differential pressure')
yyaxis right
plot(meastime,P,'o')
ylabel('tower pressure')
title('Chamber properties')

figure
plot(meastime,swin,'o')
ylabel('SWin')

figure
plot(meastime,[w h2o],'o')
ylabel('H_2O mixing ratio')
legend({'chamber','tower'})


%% main flux calculation loop

for i = 1:height(m)
    inds = timerange(m.Time(i)+minutes(1.25),m.endtime(i)-minutes(0.5)); %an index of rows that are between the start (inclusive) and end (exclusive) times
    %     figure(1)
    %     hold on
    %     plot(tt.Time(inds),tt.CH4_dry(inds),'ro') %update the plot with red dots for measurement times
    %     hold off
    dumtime = datenum(datevec(tt.Time(inds))); %convert meas. time to matlab serial date format (fractional days)
    x = dumtime-floor(dumtime); %subtract off the integer value of the date to get fractional day
    
    V = Vmisc + (S.*(m.offset(i)+hc)); % chamber volume [cm^3] incoporating offset into total chamber volume
    
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
    
    F = (10.*vol.*P(i).*(1-(w(i)/1000)).* dcdt) ./ (R.*S.*(Tfinal(i)+273.15)); %main flux equation
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
    
    F = (10.*vol.*P(i).*(1-(w(i)/1000)).* dcdt) ./ (R.*S.*(Tfinal(i)+273.15));
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
    
    %follow instructions here to install if you want to write all figures
    %to PDF https://github.com/altmany/export_fig
    %export_fig('Annandale-06-14-2022v2.pdf','-append','-painters')
    
end

%% flux figures
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

%% construct and write output file

output = m;
output.co2f = fluxes.co2;
output.ch4f = fluxes.ch4;
output.P = P;
output.RH = RHfinal;
output.Tcham = Tfinal;
output.diffpress = diffpress;
output.h2o = h2o;
output.plantflag = double(isplant);
output.Rhtower = rh;
output.SWin = swin;
output.Tairtower = ta;

timestamp = datestr(now,30); %unique timestamp at time of writing to prevent overwriting output file
fname = ['.\output\output-2022-06-14-Annandale-',timestamp,'.csv'];
mkdir('output')
writetimetable(output,fname)

%%
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




