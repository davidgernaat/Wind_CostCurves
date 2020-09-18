%% Solar_output_maps

%% Geographical potential
% CSP
CSP_GeoPot_Cellw = TechPotCellTheo{13};
CSP_GeoPot_Cellw(find(CSP_GeoPot_Cellw==0))=-99;

file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\Maps\\CSPGeoPot.asc')); 
txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.7f\nNODATA_value\t%d',nc,nr,-180,-90,0.5,-99);
dlmwrite(file,txt,'');
dlmwrite(file,CSP_GeoPot_Cellw,'-append','delimiter','\t');

% PV
PV_GeoPot_Cellw = PV_GeoPot_Cell{13};
PV_GeoPot_Cellw(find(PV_GeoPot_Cellw==0))=-99;

file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\Maps\\PVGeoPot.asc')); 
txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.7f\nNODATA_value\t%d',nc,nr,-180,-90,0.5,-99);
dlmwrite(file,txt,'');
dlmwrite(file,PV_GeoPot_Cellw,'-append','delimiter','\t');

% PVres
PVres_GeoPotw = PVres_GeoPot{13};
PVres_GeoPotw(find(PVres_GeoPotw==0))=-99;

file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\Maps\\PVresGeoPot.asc')); 
txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.7f\nNODATA_value\t%d',nc,nr,-180,-90,0.5,-99);
dlmwrite(file,txt,'');
dlmwrite(file,PVres_GeoPotw,'-append','delimiter','\t');

%% Technical potential
% CSP
AnnTechPotCellCSPw = AnnTechPotCellCSP;
AnnTechPotCellCSPw(find(AnnTechPotCellCSPw==0))=-99;

file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\Maps\\CSPTechPot.asc')); 
txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.7f\nNODATA_value\t%d',nc,nr,-180,-90,0.5,-99);
dlmwrite(file,txt,'');
dlmwrite(file,AnnTechPotCellCSPw,'-append','delimiter','\t');

% PV
AnnTechPotCellPVw = AnnTechPotCellPV;
AnnTechPotCellPVw(find(AnnTechPotCellPVw==0))=-99;

file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\Maps\\PVTechPot.asc')); 
txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.7f\nNODATA_value\t%d',nc,nr,-180,-90,0.5,-99);
dlmwrite(file,txt,'');
dlmwrite(file,AnnTechPotCellPVw,'-append','delimiter','\t');

% PVres
AnnTechPotCellPVresw = AnnTechPotCellPVres;
AnnTechPotCellPVresw(find(AnnTechPotCellPVresw==0))=-99;

file = fullfile(pathname, sprintf('\\CSP_PV_PVres\\output\\Maps\\PVresTechPot.asc')); 
txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.7f\nNODATA_value\t%d',nc,nr,-180,-90,0.5,-99);
dlmwrite(file,txt,'');
dlmwrite(file,AnnTechPotCellPVresw,'-append','delimiter','\t');

%% Economic potential

