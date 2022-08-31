function [sd] = csitime2date(year,doy,hhmm)
%CSITIME2DATE converts year, DOY, HHMM to MATLAB serial date
%   Detailed explanation goes here

hr = floor(hhmm./100);
mn = hhmm-hr*100;

m = repmat(1,size(doy));%this hack appears to work.  sent month to 1 and then use DOY for D
s = repmat(0,size(doy));
sd = datenum(year,m,doy,hr,mn,s);

end

