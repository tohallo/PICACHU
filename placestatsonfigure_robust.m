function [output] = placestatsonfigure_robust(x, y, flag,color,fid)
%placestatsonfigure_robust Robust (outlier rejection) linear regression
%   A series of fitting stats are output and plot is made
%## inputs ##
%x - x variable
%y - y variable
%flag = if flag==1, add regresion model to current figure, do nothing otherwise
%color - line color for model line
%fid - figure id to try to add regression to a figure other than the
%current one (not sure this works)

dum = [x y];
d = cleannan(dum); %clean nans
x = d(:,1); %reassign cleaned
y = d(:,2); %reassign cleaned

b = sortrows([x y],1); %sort by ascending x variable
x = b(:,1);
y = b(:,2);

% figure
% scatter(x,y)
N = length(x);
[rr,p] = corrcoef(x, y);
rs = rr(1,2)^2;
pval = p(1,2);
rootmean = rmse(x,y);

mdl = fitlm(x,y,'linear','RobustOpts','off') %built in matlab fit linear model
ci = coefCI(mdl);
ypred = predict(mdl,x);


%fprintf(fid,'%f,%f,%f,%f,%f,%f,%f\n',serdays(i),...
%output = [mdl.Coefficients.Estimate(1) ci(1,1) ci(1,2) mdl.Coefficients.Estimate(2) ci(2,1) ci(2,2)];
%pack stats structure array with linear model for output

stats.intercept =    mdl.Coefficients.Estimate(1);
stats.interceptCIlower = ci(1,1);
stats.interceptCIupper = ci(1,2);
stats.interceptPval = mdl.Coefficients.pValue(1);
stats.interceptSE = mdl.Coefficients.SE(1);
stats.slope =    mdl.Coefficients.Estimate(2);
stats.slopeCIlower = ci(2,1);
stats.slopeCIupper = ci(2,2);
stats.slopePval = mdl.Coefficients.pValue(2);
stats.slopeSE = mdl.Coefficients.SE(2);
stats.R2 = mdl.Rsquared.Ordinary;
stats.model = mdl;
output = stats;

%% add regression statistics to figure
ylims = ylim;
xlims = xlim;
text(0.1*(xlims(2)-xlims(1))+xlims(1),0.95*(ylims(2)-ylims(1))+ylims(1),sprintf('R^2=%3.2f',mdl.Rsquared.Ordinary));
text(0.1*(xlims(2)-xlims(1))+xlims(1),0.90*(ylims(2)-ylims(1))+ylims(1),sprintf('pval=%0.3g',stats.slopePval));
text(0.1*(xlims(2)-xlims(1))+xlims(1),0.85*(ylims(2)-ylims(1))+ylims(1),sprintf('slope=%0.3g',stats.slope));
text(0.1*(xlims(2)-xlims(1))+xlims(1),0.80*(ylims(2)-ylims(1))+ylims(1),sprintf('intercept=%0.3g',stats.intercept));
text(0.1*(xlims(2)-xlims(1))+xlims(1),0.75*(ylims(2)-ylims(1))+ylims(1),sprintf('RMSE=%0.3g',mdl.RMSE));

if (nargin < 4), color = 'k'; end %if no color, make it k (black)

%add model line to figure
if flag==1
    hold on
    figure(fid)
    plot((x),ypred,color,'linewidth',1.5)
    hold off
end


end

