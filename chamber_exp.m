function [C] = chamber_exp(beta,X,c0,t0)
%CHAMBER_EXP gas concentration C as an exponential function of time t
% ##inputs##
%beta - initial guess of parameters to fit
%X - x variable (time stamp)
%c0 - starting gas concentration from linear fit of deadband time
%t0 - time associated with c0 (if needs to be subtracted)
% ##output##
% C = gas concentration

%c0 is starting gas concentration as y intercept of regression of first
%10sec values (do this first in main program and send here)
%t0 is time of c0

%parameters
a = beta(1); % exponential parameter
cx = beta(2); %asymptote parameter
%obs to fit
t = X(:,1); %time

%model
%equation 1-18
%https://www.licor.com/env/support/LI-8100A/topics/deriving-the-flux-equation
C = cx + ((c0 - cx) .* exp(-a.*(t-t0)));

end

