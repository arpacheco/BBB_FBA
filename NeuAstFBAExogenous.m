%% Neuron-Astrocyte Flux Balance Analysis
% Performs flux balance analysis on the GABA-ergic model 'iNL403_GABA.mat'.
% Modifies neuron-astrocyte genome-scale metabolic model to reflect
% experimental conditions in Brain-on-a-Chip device, and tests dependence
% of neuronal GABA production on lactate and pyruvate from the environment,
% and on glucose from the environment. Metabolic model sourced from: Lewis
% N.E. et al., Nat. Biotech. 2010. FBA optimization performed using COBRA
% Toolbox v2.0, Schellenberger et al. Nat. Protocols, 2011. Solver: Gurobi
% 6.
%
% Alan R. Pacheco, Boston University Graduate Program in Bioinformatics
% 17-Jul 2016

%% Initialize the COBRA Toolbox
% initCobraToolbox

%% Load model
filename = 'iNL403_GABA.mat';
S = load(filename,'-mat');
Name = whos('-file',filename);
model = Name.name;
model = S.(model);

%% Run FBA to get ATP demand fluxes
FBAsoln = optimizeCbModel(model,'max',0);
ATP_As_min = FBAsoln.x(10); % astrocyte
ATP_Ne_min = FBAsoln.x(603); % neuron

%% Correct gln-glu-gaba exchange directionalities
% correct directionality of gaba exchange
model.S(find(model.S(:,668)),668) = -model.S(find(model.S(:,668)),668);
model.lb(668)=0;
model.ub(668)=1000;
model.lb(536)=0;

%% Close interstitial lac and pyruvate transport and make lac and pyr come from extracellular to neuron
% constrain lac and pyr to neuron from interstitial
model.lb([880 983])=0;
model.ub([880 983])=0;

% constrain lac and pyr from astrocyte to interstitial
model.lb([567 580])=0;
model.ub([567 580])=0;

% Create new extracellular -> neuron lactate reaction
model.S(872,1068) = 1; % lac-L [cN]
model.S(561,1068) = 1; % h[cN]
model.S(168,1068) = -1; % lac-L[e]
model.S(161,1068) = -1; % h[e]
model.rxns{1068} = 'L-LACt2r_neu';
model.rxnNames{1068} = 'L-lactate reversible transport via proton symport Neuron';
model.c(1068) = 0;
model.rev(1068) = 1;
model.lb(1068) = -1000;
model.ub(1068) = 1000;

% Create new extracellular -> neuron pyruvate reaction
model.S(708,1069) = 1; % pyr [cN]
model.S(561,1069) = 1; % h[cN]
model.S(181,1069) = -1; % pyr[e]
model.S(161,1069) = -1; % h[e]
model.rxns{1069} = 'PYRt2r_neu';
model.rxnNames{1069} = 'pyruvate reversible transport via proton symport Neuron';
model.c(1069) = 0;
model.rev(1069) = 1;
model.lb(1069) = -1000;
model.ub(1069) = 1000;

%% Change objective to be GABA output by neuron but still keep ATP demand minimum
model.c(:) = 0;
obj = [823];
model.c(obj)=1;
model.lb(10)=ATP_As_min; % Astrocyte ATP demand
model.lb(603)=ATP_Ne_min; % Neuron ATP demand

%% Vary global intake of glc, lac, pyr, and record effect on GABA production
lac_lb_orig = model.lb(45);
pyr_lb_orig = model.lb(58);
glc_lb_orig = model.lb(33);
met_up_range = [-0.6:0.001:0]';
model.ub(45)=0;
model.ub(58)=0;
model.ub(33)=0;

% vary glucose with constant lactate and pyruvate
for i = 1:length(met_up_range)
    model.lb(33) = met_up_range(i);
    FBAsoln = optimizeCbModel(model,'max',0);
    if numel(FBAsoln.x)==0
        gaba_exc_glc_1(i) = 0;
    else
        gaba_exc_glc_1(i) = FBAsoln.x(823);
    end
end
model.lb(33) = glc_lb_orig;

% vary glucose without lactate or pyruvate
model.lb(45)=0;
model.lb(58)=0;
for i = 1:length(met_up_range)
    model.lb(33) = met_up_range(i);
    FBAsoln = optimizeCbModel(model,'max',0);
    if numel(FBAsoln.x)==0
        gaba_exc_glc_2(i) = 0;
    else
        gaba_exc_glc_2(i) = FBAsoln.x(823);
    end
end
model.lb(33) = glc_lb_orig;

% vary lactate and pyruvate together with constant glucose
for i = 1:length(met_up_range)
    model.lb(45) = met_up_range(i);
    model.lb(58) = met_up_range(i);
    FBAsoln = optimizeCbModel(model,'max',0);
    if numel(FBAsoln.x)==0
        gaba_exc_lacpyr_1(i) = 0;
    else
        gaba_exc_lacpyr_1(i) = FBAsoln.x(823);
    end
end
model.lb(45) = lac_lb_orig;
model.lb(58) = pyr_lb_orig;

% vary lactate and pyruvate together without glucose
model.lb(33) = 0;
for i = 1:length(met_up_range)
    model.lb(45) = met_up_range(i);
    model.lb(58) = met_up_range(i);
    FBAsoln = optimizeCbModel(model,'max',0);
    if numel(FBAsoln.x)==0
        gaba_exc_lacpyr_2(i) = 0;
    else
        gaba_exc_lacpyr_2(i) = FBAsoln.x(823);
    end
end

dataMat = [gaba_exc_glc_1',gaba_exc_glc_2',gaba_exc_lacpyr_1',gaba_exc_lacpyr_2'];

%% Plotting
close all

% Plot GABA as a function of overall glc,lac,pyr intake
figure
plot(-met_up_range,gaba_exc_glc_1,'Linewidth',4)
hold on
plot(-met_up_range,gaba_exc_glc_2,'Linewidth',2)
hold on
plot(-met_up_range,gaba_exc_lacpyr_1,'Linewidth',2)
hold on
plot(-met_up_range,gaba_exc_lacpyr_2,'Linewidth',2)
title('Production of GABA')
legend('Varied glucose uptake (with constant lactate and pyruvate)','Varied glucose uptake (no lactate or pyruvate)','Varied lactate and pyruvate uptake (with constant glucose)','Varied lactate and pyruvate uptake (no glucose)','Location','southeast')
xlabel('Metabolite uptake flux (umol/gWB/min)')
ylabel('Synaptic cleft GABA flux (umol/gWB/min)')