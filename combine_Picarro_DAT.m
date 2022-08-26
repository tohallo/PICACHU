function [T2] = combine_Picarro_DAT(currentfolder,destinationfolder,starttime,endtime)
%COMBINE_Picarro_DAT Combines Picarro .dat files into one table
%inputs
%currentfolder - path of referring function (pwd)
%destinationfolder- directory relative to where main .m file lives
%starttime - starting time to look for files (matlab serial date)
%endtime - ending time to look for files (matlab serial date)
%output
%T2 - concatenated table in matlab timetable format

cd(destinationfolder) % cd to location to write concatenated table
giant = table(); %initialize empty table
d = dir('*.dat'); %list of all dat files
for i = 1:length(d) %loop through dat files
    dum=readtable([d(i).name]); % read file as table into variable dum
    giant = tblvertcat(giant,dum); % a function from the matlab file exchange that concatenates data tables, removing headers (adds dum to bottom of giant)
    disp(i) %print the i value of the loop
end

%DT1 = datetime(datevec(ameriflux2serdate(giant.TIMESTAMP_START)));
%DT1=datetime(join([string(giant.date) string(giant.time)]));
DT1=datetime(join([string(giant.DATE) string(giant.TIME)])); %create date vector in matlab datetime format
% giant.datetime = csitime2date(giant.YEAR,giant.JULIEN,giant.Time);
T = table2timetable(giant,'RowTimes',DT1); %convert to matlab timetable format
%T2 = retime(sortrows(T),'regular','nearest','TimeStep',minutes(30));
dum = sort(T.Properties.RowTimes); %sort in case of duplicates
%newtimevector = floor(datenum(dum(1))):1/48:(floor(datenum(dum(end)))+1)-1/48; %round to midnight of starting and ending days
newtimevector = starttime:1/(1440*60):endtime; %make new time vector (one second time step) using input arguments in matlab datenum format (not currently used)

% T2 =retime(sortrows(T),datetime(datevec(newtimevector))','fillwithmissing'); %interpolate tabel to new time stamp (not using because do not want to use interpolated values
% T2(T2==-7999)=NaN;
T2 = T; %just use T for now, this is output variable (combined table)  
cd(currentfolder) %return to referrng directory
end




