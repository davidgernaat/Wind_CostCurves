%% Solar Theo+Tech+Econ potential per country


ISOCountryNumbers = ISOGDP(:,1);
ISOCountryNumbers(end+1) = 0;
ISOCountryNames = CountryNames';
ISOCountryNames{end+1} = 'World';

conv=1e-9; % kWh --> TWh

TheoPotWon = round(CTechPotWindTheo{13}*conv)';
TheoPotWoff = round(CTechPotWindTheoo{13}*conv)';

GeoPotWon = round(CTechPotWindGeo{13}*conv)';
GeoPotWoff = round(CTechPotWindGeoo{13}*conv)';

TechPotWon = round(CTechPotWind{13}*conv)';
TechPotWoff = round(CTechPotWindo{13}*conv)';

EconPotWon = round(CTechPotWindEcon{13}*conv)';
EconPotWoff = round(CTechPotWindEcono{13}*conv)';

T = table(ISOCountryNumbers, ISOCountryNames,...
    TheoPotWon,TheoPotWoff,...
    GeoPotWon,GeoPotWoff,...
    TechPotWon,TechPotWoff,...
    EconPotWon,EconPotWoff); 

file = fullfile(root, sprintf('\\output\\Wind_Potentials_per_Country.csv'));

writetable(T,file,'Delimiter',',')