function [CostCurveSmthOnshore, CostCurveSmthOffshore, ...
    MaxProdOnshore, MaxProdOffshore, ...
    LFCurveSmthOnshore, LFCurveSmthOffshore, ...
    TechPotCellTheo, TechPotCellGeo, TechPotCell, COE] = Wind_cc_ISIMIP(root, V)

% V=CMc;
%Windfile is set to 1 so the NASA wind files do not overwrite ISIMIP V

%% CSP PV and PVresidential cost supply curve generator
% clear all
% root = 'Y:\ontwapps\Timer\Users\David\Pojects\Offshore wind\Cost curve prog\Wind_CC';

fname = sprintf('%s\\input\\input_data_onshore_ISIMIP.mat', root);
load(fname)
fname = sprintf('%s\\input\\input_data_offshore_ISIMIP.mat', root);
load(fname)

[nr nc] = size(IRegion);
%% Settings
% 0=no/1=yes
setoutput=0;                % Do you want output?
outputcountrypots=0;        % Do you want to output country potentials?
setOnshorePL=1;             % Do you want to include powerline costs?
setEconCeil=0;              % Do you want to add a cost ceiling for economic potential?
hubheightcorr=0;            % do you want the NASA wind data corrected to 10m hub height? Better comparability with hoogwijk

glctfile=2;                 % 1=old file from hoogwijk / 2=newer glct file (used in Solar_cc)
windfile=1;                 % 1=Hoogwijk winds (m/s) / 2=NASA wind speeds (m/s)
areafile=2;                 % 1=Hoogwijk file (much smaller in higher lats) / 2=self calculated
vcorrection=0;              % 0=no V corrections / 1=V correction similar to Hoogwijk >10m/s=8m/s

CostCeil=0.15;               % Level of cost ceiling $ / kWh

cmap = jet(256);
cmap(1,:) = [1 1 1];

ir=13; %first imagemask row number compatible with MyM
ic=284; %first imagemask column number compatible with MyM

%% Literature numbers
fprintf('\nGlobal Technical Potential Onshore Hoogwijk %.2f  PWh', 96);
fprintf('\nGlobal   Technical Potential Onshore Seurek %.2f PWh', 557);
fprintf('\nGlobal       Technical Potential Onshore Lu %.2f PWh', 690);
fprintf('\nGlobal    Technical Potential Offshore NREL %.2f PWh', 276);
fprintf('\nGlobal  Technical Potential Offshore Seurek %.2f PWh\n', 315);

%% Prep
disp('Prep')
Wind_prep

%% Technical potential 1
disp('Technical potential 1')
Wind_TechPot1
% Calculation of technical potential parameters
% We first correct for the height using the roughness length
% The roughness length depends on the land use type
% furthermore, the wind speed is extrapolated to the wind speed at hubheight.
% The hubheight is exogenous determined and increases over time.

%% Geographical potential
disp('Geographical potential')
Wind_GeoPot
% The site potential is determined by land cover data and population density data.

%% Technical potential 2
disp('Technical potential 2')
Wind_TechPot2
% Determines power output per cell
%%
%Onshore
fprintf('\nGlobal              Technical Potential Onshore %.2f PWh', RegTechPotWind{13}(27)*1e-12);
fprintf('\nGlobal       No cut Technical Potential Onshore %.2f PWh', RegTechPotWindNoCut{13}(27)*1e-12);
fprintf('\nGlobal  Theoretical Technical Potential Onshore %.2f PWh', RegTechPotWindTheo{13}(27)*1e-12);

%Offshore
fprintf('\n\nGlobal             Technical Potential Offshore %.2f PWh', RegTechPotWindo{13}(27)*1e-12);
fprintf('\nGlobal      No cut Technical Potential Offshore %.2f PWh', RegTechPotWindNoCuto{13}(27)*1e-12);
fprintf('\nGlobal Theoretical Technical Potential Offshore %.2f PWh\n', RegTechPotWindTheoo{13}(27)*1e-12);

%% Economic potential
disp('Economic potential')
Wind_EconPot
% Determines economic potential
%%
if setEconCeil
    %Onshore
    fprintf('\nGlobal Economic <0.5 $/kWh Technical Potential Onshore %.2f PWh', RegTechPotWindEcon{13}(27)*1e-12);
    %Offshore
    fprintf('\nGlobal Economic <0.5 $/kWh Technical Potential Offshore %.2f PWh\n', RegTechPotWindEcono{13}(27)*1e-12);
end

%% Costsupply
disp('Cost Supply')
Wind_CostSupply

%% Output country potentials
if outputcountrypots
    disp('Country potentials output')
    Wind_PotsPerCountry
end
