
clear model data params options
load Malariadata-pv.mat

%%
model.ssfun = @Malaria2;

%%
% All parameters are constrained to be positive. The initial
% concentrations are also unknown and are treated as extra parameters.
params = {
    {'ei1', 13,  10,20, 13, 1}  %Time from exposed to infected
    {'i1i2',  115, 90,200,100,10}%Time from symptomatic to asymptomatic
    {'s2s1', 10,   0,  20, 10, 1 }  %Time of recovery from subpatent infection to susceptible
    {'sf', 0.5,   0,  1 }  %Infectivity of subpatent infected people
    {'ps', 0.08,   0.06,  0.12} %reporting rate 
    {'si2', 0.5,   0,  1 }  %Superinfection fraction, from S2 to I1
    {'ts', 0.9,   0.7,  1 }  %Treatment success 
    {'s1', 0.5, 0,1}
    {'s2', 0.5, 0,1}
    {'s3', 0.5, 0,1}
    {'s4', 0.5, 0,1}
    {'s5', 0.5, 0,1}
    {'s6', 0.5, 0,1}
    {'iq1',0.7,0,1}
    {'qi1',210,0,1000,210,30}
    {'betanino3',0.01,0,Inf}
    {'spray',0.001,0,Inf}
    {'net',0.01,0,Inf}
    {'HI1_0', 8E-4, 0, 0.1}
    {'HS2_0', 5E-2, 0, 0.1}  
    };

%%
model.S20 = [2.5];
model.N0  = [0];

%%
% First generate an initial chain.
options.nsimu = 1000000;
options.stats = 1;
[results, chain, s2chain]= mcmcrun(model,data,params,options);
%%
options.nsimu = 1000000;
options.stats = 1;
[results, chain, s2chain] = mcmcrun(model,data,params,options, results);

%%
figure
mcmcplot(chain,[],results); %,'pairs'
figure
mcmcplot(chain,[],results,'denspanel',2);

%%
chainstats(chain,results)

%%
% In order to use the |mcmcpred| function we need
% function |modelfun| with input arguments given as
% |modelfun(xdata,theta)|. We construct this as an anonymous function.

modelfun = @(d,th) Malaria3(d(:,1),th,th(end-1:end),d);

nsample = 10000;
results.sstype = 1;
out = mcmcpred(results,chain,s2chain,data.xdata,modelfun,nsample);
figure
mcmcpredplot(out);

for i = 1:2
 hold on
 subplot(3,1,i)
 hold on
 plot(data.ydata(:,1),data.ydata(:,i+1),'s');
 hold off
end

fit=out.predlims{1,1}{1,1}(2,:)';
r2=1-sum((data.ydata(:,2)-fit).^2)/sum((data.ydata(:,2)-mean(data.ydata(:,2))).^2)
 