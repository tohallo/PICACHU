function [Xclean] = cleannan(X)
X2 = X;

X2(any(isnan(X2),2),:) = []; %clean out nans

Xclean = X2;
