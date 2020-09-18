%% Run ISIMIPs
clear all

root = 'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER\Wind_CC';
root2 = 'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER';
root_data = 'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER\Wind_CC\input\ISIMIP\Wind';

% 1. tas files : are in Kelvin ...so 273.15 must be subtracted from the values to get degree Celsius (K-  273.15=oC)
% 2. sfWind files: are in m s-1
% 3. rsds files : are in W m-2

fname = sprintf('%s\\input\\input_data_offshore.mat', root);
load(fname)

%% File list
disp('Read file list')
fnames = dir(root_data);

c=0;
for j=1:numel(fnames)
    if j<3; continue; end;
    c=c+1;
    pathnames{c} = sprintf('%s\\%s',root_data,(fnames(j).name));
end

c=0;
% clear names
% clc
for j=3:numel(fnames)
    c=c+1;
    
    if c<6
        %fprintf('%s\n',fnames(j).name(13:70))
        names{c} = fnames(j).name(13:70);
    elseif c>5 && c<11
        %fprintf('%s\n',fnames(j).name(13:70))
        names{c} = fnames(j).name(13:70);
    elseif c>10 && c<15
        %fprintf('%s\n',fnames(j).name(13:70))
        names{c} = fnames(j).name(13:70);
    elseif c>14
        %fprintf('%s\n',fnames(j).name(13:65))
        names{c} = fnames(j).name(13:65);
    end
    
    names_split{c} =strsplit(names{c},'_');
    
    %fprintf('%s %s %s\n',names_split{c}{1},names_split{c}{2},names_split{c}{5})
end

for i=1:numel(names_split)
    GCMID{i}=names_split{i}{1};
    RCPID{i}=names_split{i}{2};
    TIMEID{i}=names_split{i}{5};
end

%% step 1: read netcdf

i=6; %Hadgem historical
% i=9; %Hadgem RCP60

for i=1:numel(names)
    fprintf('Reading file #%d of %d\n',i, numel(names))
    clear CM CMc
    
    %% Read netcdf
    % finfo = ncinfo(names{i});
    % finfo.Variables.Name
    
    CM_data = ncread(pathnames{i},'sfcWind');
    
    %% step 2: convert to correct format
    disp('Re-formatting')
    
    for m=1:12
        for r=1:360
            for c=1:720
                CM{m}(r,c) = CM_data(c,r,m);
            end
        end
    end
    
    for r=1:360
        for c=1:720
            for m=1:12
                dt(m) = CM{m}(r,c);
            end
            CM{13}(r,c) = mean(dt);
        end
    end
    
    %% step 3: converting units
    disp('Converting')
    
    for m=1:13
        for r=1:360
            for c=1:720
                CMc{m}(r,c) = CM{m}(r,c) * 1; % data is already in m/s
            end
        end
    end
    
    % figure(1);clf;imagesc(CMc{13});colorbar
    
    %% step 4: running Solar_cc with Climate model data
    disp('Running Cost-supply')
    
    [CostCurveSmthOnshore{i}, CostCurveSmthOffshore{i}, ...
        MaxProdOnshore{i}, MaxProdOffshore{i}, ...
        LFCurveSmthOnshore{i}, LFCurveSmthOffshore{i}, ...
        TechPotCellTheo{i}, TechPotCellGeo{i}, TechPotCell{i}, COE{i}] = Wind_cc_ISIMIP(root,CMc);
    
    %% step 5: Deal with IAM output
    
    % POLES
    C2R_fname = sprintf('%s\\input\\ISIMIP\\Modelregionallocation\\Country_to_POLESregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\POLES\\Wind\\Wind_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,TechPotCell{i}{13},COE{i}{13},GCMID{i},RCPID{i},TIMEID{i})
    modelregwindoutput(root,C2R_fname,scenlib_IAM,CMc,GCMID{i},RCPID{i},TIMEID{i})

    % REMIND
    C2R_fname = sprintf('%s\\input\\ISIMIP\\Modelregionallocation\\Country_to_REMINDregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\REMIND\\Wind\\Wind_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,TechPotCell{i}{13},COE{i}{13},GCMID{i},RCPID{i},TIMEID{i})
    
    % TIAM
    C2R_fname = sprintf('%s\\input\\ISIMIP\\Modelregionallocation\\Country_to_TIAMregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\TIAM\\Wind\\Wind_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,TechPotCell{i}{13},COE{i}{13},GCMID{i},RCPID{i},TIMEID{i})
    
    % MESSAGE
    C2R_fname = sprintf('%s\\input\\ISIMIP\\Modelregionallocation\\Country_to_MESSAGEregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\MESSAGE\\Wind\\Wind_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,TechPotCell{i}{13},COE{i}{13},GCMID{i},RCPID{i},TIMEID{i})

    % COFFEE
    C2R_fname = sprintf('%s\\input\\ISIMIP\\Modelregionallocation\\Country_to_COFFEEregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\COFFEE\\Wind\\Wind_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,TechPotCell{i}{13},COE{i}{13},GCMID{i},RCPID{i},TIMEID{i})
    
    % TIAM31
    C2R_fname = sprintf('%s\\input\\ISIMIP\\Modelregionallocation\\Country_to_TIAM31Rregion_IsoCode.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\TIAM31R\\Wind\\Wind_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,TechPotCell{i}{13},COE{i}{13},GCMID{i},RCPID{i},TIMEID{i})
    
    % TIAM country division
    C2R_fname = sprintf('%s\\input\\ISIMIP\\Modelregionallocation\\GISO_country_division.csv',root);
    scenlib_IAM = sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\TIAM_country\\Wind\\Wind_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
    
    modelregoutput(root,C2R_fname,scenlib_IAM,TechPotCell{i}{13},COE{i}{13},GCMID{i},RCPID{i},TIMEID{i})
    
    
    %% step 6: Deal with output
    disp('Output')
    
    scenlib = sprintf('Y:\\Kennisbasis\\IMAGE\\model\\users\\david\\Pojects\\ISIMIP2E\\ISIMIP2E\\2_TIMER\\scenlib\\TIMER_2015\\ISIMIP2E\\maps\\Wind\\Wind_%s_%s_%s',GCMID{i},RCPID{i},TIMEID{i});
    
    matpath = fullfile(scenlib, sprintf(''));
    if ~isdir(matpath)
        mkdir(matpath);
    end
    pathname = fileparts(scenlib);
    
    % CostCurves
    
    % Onshore
%     file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\curves\\CostCurveSmthWON.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.WONCostCurveSmth[27](oo) = [');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,CostCurveSmthOnshore{i},'-append','delimiter','\t');
%     dlmwrite(file,']','-append');
%     dlmwrite(file,';','-append');
%     
%     % Offshore
%     file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\curves\\CostCurveSmthWOFF.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.WOFFCostCurveSmth[27](oo) = [');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,CostCurveSmthOffshore{i},'-append','delimiter','\t');
%     dlmwrite(file,']','-append');
%     dlmwrite(file,';','-append');
%     
%     % MaxProd
%     
%     % Onshore
%     file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\curves\\MaxProdWON.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.MaxProdOnshore[27] = ');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,MaxProdOnshore{i},'-append','delimiter','\t');
%     dlmwrite(file,';','-append');
%     
%     % Offshore
%     file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\curves\\MaxProdWOFF.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.MaxProdOffshore[27] = ');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,MaxProdOffshore{i},'-append','delimiter','\t');
%     dlmwrite(file,';','-append');
%     
%     % Load factor decline curves
%     
%     % Onshore
%     file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\curves\\LoadCurveSmthWON.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.LFCurveSmthOnshore[27](oo) = [');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,LFCurveSmthOnshore{i},'-append','delimiter','\t');
%     dlmwrite(file,']','-append');
%     dlmwrite(file,';','-append');
%     
%     % Offshore
%     file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\curves\\LoadCurveSmthWOFF.dat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt=sprintf('real main.LFCurveSmthOffshore[27](oo) = [');
%     dlmwrite(file,txt,'');
%     dlmwrite(file,LFCurveSmthOffshore{i},'-append','delimiter','\t');
%     dlmwrite(file,']','-append');
%     dlmwrite(file,';','-append');
    
    % Maps
    txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.2f\nNODATA_value\t%d',720,360,-180,-90,0.5,-99);
    
    % Wind theoretical
    file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\TheoPotWind.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    dlmwrite(file,TechPotCellTheo{i}{13},'-append','delimiter',' '); % kWh / cell / y
    %figure(1);clf;imagesc(TechPotCellTheo{i}{13});axis image;colorbar
    
    % Wind geographic
    file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\GeoPotWind.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    TechPotCellGeo{i}{13}(OffDis>5)=-99;
    dlmwrite(file,TechPotCellGeo{i}{13},'-append','delimiter',' '); % kWh / cell / y
    %figure(1);clf;imagesc(TechPotCellGeo{i}{13});axis image;colorbar
    
    % Wind technical
    file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\TechPotWind.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    TechPotCell{i}{13}(OffDis>5)=-99;
    dlmwrite(file,TechPotCell{i}{13},'-append','delimiter',' '); % kWh / cell / y
    %figure(2);clf;imagesc(TechPotCell{i}{13});axis image;colorbar
    
%     if i==6 %for show PhD mid chapter hadgem techpot
%         % Wind technical
%         file = fullfile(pathname, sprintf('\\xPot_Mean_change_maps\\Techpot_PhD_midchapter\\data\\Wind_%s_%s_%s_TechPotWind.asc',GCMID{i},RCPID{i},TIMEID{i}));
%         dt=TechPotCell{i}{13}*3.6e-9;
%         dt(OffDis>5)=-99;
%         fprintf('Wind techpot %0.2f EJ\n',sum(dt(dt~=-99))*1e-3)
%         dt(dt>150)=150;
%         dlmwrite(file,txt,'');
%         dlmwrite(file,dt,'-append','delimiter',' '); % PJ / cell / y
%         %figure(2);clf;imagesc(dt);axis image;colorbar
%         
%         % Wind economic
%         file = fullfile(pathname, sprintf('\\xPot_Mean_change_maps\\Techpot_PhD_midchapter\\data\\Wind_%s_%s_%s_EconPotWind.asc',GCMID{i},RCPID{i},TIMEID{i}));
%         dt=TechPotCell{i}{13}*3.6e-9;
%         dt(OffDis>5)=-99;
%         dt(COE{i}{13}>0.1)=0;
%         fprintf('Wind econpot %0.2f EJ\n',sum(dt(dt~=-99))*1e-3)
%         dt(dt>150)=150;
%         dlmwrite(file,txt,'');
%         dlmwrite(file,dt,'-append','delimiter',' '); % PJ / cell / y
%         %figure(2);clf;imagesc(dt);axis image;colorbar
%     end
    
    %%
    % Wind COE
    file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\COEwind.asc',GCMID{i},RCPID{i},TIMEID{i}));
    dlmwrite(file,txt,'');
    COE{i}{13}(OffDis>5)=-99;
    COE{i}{13}(find(isnan(COE{i}{13}(:))))=-99;
    COE{i}{13}(find(isinf(COE{i}{13}(:))))=-99;
    dlmwrite(file,COE{i}{13},'-append','delimiter',' '); % $2005/kWh
    %figure(1);clf;imagesc(COE{i}{13});axis image;colorbar
    
    %% input sce-file
    % Keep format. in the sce-file the format is correct
%     file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\scenario_input.sce',GCMID{i},RCPID{i},TIMEID{i}));
%     txt1 = sprintf('DIRECTORY("../scenlib/$1/curves");');
%     txt2 = sprintf('FILE("CostCurveSmthWON.dat","r")	=  main.em.ep.isod.CostCurveSmthWON;');
%     txt3 = sprintf('FILE("CostCurveSmthWOFF.dat","r")	=  main.em.ep.isod.CostCurveSmthWOFF;');
%     txt4 = sprintf('FILE("MaxProdWON.dat","r")              =  main.em.ep.isod.MaxProdWON;');
%     txt5 = sprintf('FILE("MaxProdWOFF.dat","r")             =  main.em.ep.isod.MaxProdWOFF;');
%     txt6 = sprintf('FILE("LoadCurveSmthWON.dat","r")	=  main.em.ep.isod.LoadCurveSmthWON;');
%     txt7 = sprintf('FILE("LoadCurveSmthWOFF.dat","r")	=  main.em.ep.isod.LoadCurveSmthWOFF;');
%     
%     dlmwrite(file,txt1,'')
%     dlmwrite(file,txt2,'-append','delimiter','')
%     dlmwrite(file,txt3,'-append','delimiter','')
%     dlmwrite(file,txt4,'-append','delimiter','')
%     dlmwrite(file,txt5,'-append','delimiter','')
%     dlmwrite(file,txt6,'-append','delimiter','')
%     dlmwrite(file,txt7,'-append','delimiter','')
%     
%     %% output sce-file
%     file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\scenario_output.sce',GCMID{i},RCPID{i},TIMEID{i}));
%     txt1 = sprintf('DIRECTORY("../outputlib/$1");');
%     txt2 = sprintf('FILE("ElecProdSpec.out", "w")			= main.em.ep.opr.ElecProdSpec2;');
%     txt3 = sprintf('FILE("GCap.out","w")                            = main.em.ep.gcp.GCap;');
%     txt4 = sprintf('FILE("CO2Spec.out","w") 			= main.em.CO2Spec;');
%     txt5 = sprintf('FILE("tpesEXT.out","w")                         = main.em.TPEStotalEXT;');
%     txt6 = sprintf('FILE("TotalCostPerkWhNew.out", "w")		= main.em.ep.epc.TotalCostPerkWhNew;');
%     txt7 = sprintf('FILE("TotalCostPerkWhAvg.out", "w")		= main.em.ep.epc.TotalCostPerkWhAvg;');
%     txt8 = sprintf('FILE("StorLossISOTot.out", "w")			= main.em.ep.isp.StorLossISOTot;');
%     txt9 = sprintf('FILE("StorLossFracISO.out", "w")		= main.em.ep.isp.StorLossFracISO;');
%     txt10 = sprintf('FILE("ISODeple.out", "w")                       = main.em.ep.isod.ISODeple;');
%     txt11 = sprintf('FILE("ISOAddOnNew.out", "w")			= main.em.ep.epc.ISOAddOnNew;');
%     txt12 = sprintf('FILE("ISOAddOnAvg.out", "w")			= main.em.ep.epc.ISOAddOnAvg;');
%     txt13 = sprintf('FILE("TotalCostperkWhAvgPVres.out","w")         = main.em.ep.epc.TotalCostperkWhAvgPVres;');
%     txt14 = sprintf('FILE("BLF_Prod.out","w") 			= main.em.bio.BCap.BLFProdTot;');
%     txt15 = sprintf('FILE("BSF_Prod.out","w") 			= main.em.bio.BCap.BSFProdTot;');
%     txt16 = sprintf('FILE("tpes.out","w")                         	= main.em.TPEStotal;');
%     
%     dlmwrite(file,txt1,'')
%     dlmwrite(file,txt2,'-append','delimiter','')
%     dlmwrite(file,txt3,'-append','delimiter','')
%     dlmwrite(file,txt4,'-append','delimiter','')
%     dlmwrite(file,txt5,'-append','delimiter','')
%     dlmwrite(file,txt6,'-append','delimiter','')
%     dlmwrite(file,txt7,'-append','delimiter','')
%     dlmwrite(file,txt8,'-append','delimiter','')
%     dlmwrite(file,txt9,'-append','delimiter','')
%     dlmwrite(file,txt10,'-append','delimiter','')
%     dlmwrite(file,txt11,'-append','delimiter','')
%     dlmwrite(file,txt12,'-append','delimiter','')
%     dlmwrite(file,txt13,'-append','delimiter','')
%     dlmwrite(file,txt14,'-append','delimiter','')
%     dlmwrite(file,txt15,'-append','delimiter','')
%     dlmwrite(file,txt16,'-append','delimiter','')
%     
%     %% settings bat-file
%     file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\scenario_settings.bat',GCMID{i},RCPID{i},TIMEID{i}));
%     txt1 = sprintf('REM Scenarios settings');
%     dlmwrite(file,txt1,'')
    
end

%% Time dependent output for TIMER with climate change impacts

tv=[1971 2000 2050 2085 2100];
% for time depended curves
jv=[1 1 2 3 3 1 1 4 5 5 ... % GFDL-ESM2M RCP26 hist 2050 2085 2100; RCP60 hist 2050 2085 2100
    6 6 7 8 8 6 6 9 10 10 ... % HadGEM-ES RCP26 hist 2050 2085 2100; RCP60 hist 2050 2085 2100
    11 11 12 13 13 11 11 14 15 15 ... %IPSL-CM5A-LR RCP26 hist 2050 2085 2100; RCP60 hist 2050 2085 2100
    16 16 17 18 18 16 16 19 20 20]; % MIROC5 RCP26 hist 2050 2085 2100; RCP60 hist 2050 2085 2100

j=0;
k=0;
l=0;
%Prep order curves
for i=1:numel(jv)
    j=j+1;
    if j>5; j=1; end
    if j==1; k=k+1; l=0; end
    %fprintf('t=%d j=%d jv=%d k=%d l=%d\n',i,j,jv(i),k,l)
    
    for c=1:101 %cost-curve steps
        for R=1:27 %region
            l=l+1;
            if c==1 && R==1
                CostCurveSmthOnshoret{k}(l)=tv(j);
                CostCurveSmthOffshoret{k}(l)=tv(j);
                LFCurveSmthOnshoret{k}(l)=tv(j);
                LFCurveSmthOffshoret{k}(l)=tv(j);
                l=l+1;
                CostCurveSmthOnshoret{k}(l)=CostCurveSmthOnshore{jv(i)}(c,R+1);
                CostCurveSmthOffshoret{k}(l)=CostCurveSmthOffshore{jv(i)}(c,R+1);
                LFCurveSmthOnshoret{k}(l)=LFCurveSmthOnshore{jv(i)}(c,R+1);
                LFCurveSmthOffshoret{k}(l)=LFCurveSmthOffshore{jv(i)}(c,R+1);
            else
                CostCurveSmthOnshoret{k}(l)=CostCurveSmthOnshore{jv(i)}(c,R+1);
                CostCurveSmthOffshoret{k}(l)=CostCurveSmthOffshore{jv(i)}(c,R+1);
                LFCurveSmthOnshoret{k}(l)=LFCurveSmthOnshore{jv(i)}(c,R+1);
                LFCurveSmthOffshoret{k}(l)=LFCurveSmthOffshore{jv(i)}(c,R+1);
            end
            
        end
        
    end
end

% prep order maxprod
j=0;
k=0;
l=0;
for i=1:numel(jv)
    j=j+1;
    if j>5; j=1; end
    if j==1; k=k+1; l=0; end
    %fprintf('t=%d j=%d jv=%d k=%d l=%d\n',i,j,jv(i),k,l)
    
    for R=1:27 %region
        l=l+1;
        if R==1
            MaxProdOnshoret{k}(l)=tv(j);
            MaxProdOffshoret{k}(l)=tv(j);
            l=l+1;
            MaxProdOnshoret{k}(l)=MaxProdOnshore{jv(i)}(R);
            MaxProdOffshoret{k}(l)=MaxProdOffshore{jv(i)}(R);
        else
            MaxProdOnshoret{k}(l)=MaxProdOnshore{jv(i)}(R);
            MaxProdOffshoret{k}(l)=MaxProdOffshore{jv(i)}(R);
        end
    end
end


% Writing
GCMIDt = {'GFDL-ESM2M' 'HADGEM2-ES' 'IPSL-CM5A-LR' 'MIROC5'};
RCPIDt = {'RCP26' 'RCP60'};
c=0;
for i=1:4
    for j=1:2
        c=c+1;
        scenlib = sprintf('Y:\\Kennisbasis\\IMAGE\\model\\users\\david\\Pojects\\ISIMIP2E\\ISIMIP2E\\2_TIMER\\scenlib\\TIMER_2015\\ISIMIP2E\\Wind_%s_%s',GCMIDt{i},RCPIDt{j});
        matpath = fullfile(scenlib, sprintf('\\curves'));
        if ~isdir(matpath); mkdir(matpath); end
        pathname = fileparts(scenlib);
        
        % Onshore cost
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\CostCurveSmthWON.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.WONCostCurveSmth[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,CostCurveSmthOnshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % Offshore cost
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\CostCurveSmthWOFF.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.WOFFCostCurveSmth[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,CostCurveSmthOffshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % Onshore MaxProd
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\MaxProdWON.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.MaxProdOnshore[27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,MaxProdOnshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % Offshore MaxProd
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\MaxProdWOFF.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.MaxProdOffshore[27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,MaxProdOffshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % Onshore load
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\LoadCurveSmthWON.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.LFCurveSmthOnshore[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,LFCurveSmthOnshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % Offshore load
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\LoadCurveSmthWOFF.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.LFCurveSmthOffshore[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,LFCurveSmthOffshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        %% input sce-file
        % Keep format. in the sce-file the format is correct
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\scenario_input.sce',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('DIRECTORY("../scenlib/$1/curves");');
        txt2 = sprintf('FILE("CostCurveSmthWON.dat","r")	=  main.em.ep.isod.CostCurveSmthWON;');
        txt3 = sprintf('FILE("CostCurveSmthWOFF.dat","r")	=  main.em.ep.isod.CostCurveSmthWOFF;');
        txt4 = sprintf('FILE("MaxProdWON.dat","r")              =  main.em.ep.isod.MaxProdWON;');
        txt5 = sprintf('FILE("MaxProdWOFF.dat","r")             =  main.em.ep.isod.MaxProdWOFF;');
        txt6 = sprintf('FILE("LoadCurveSmthWON.dat","r")	=  main.em.ep.isod.LoadCurveSmthWON;');
        txt7 = sprintf('FILE("LoadCurveSmthWOFF.dat","r")	=  main.em.ep.isod.LoadCurveSmthWOFF;');
        
        dlmwrite(file,txt1,'')
        dlmwrite(file,txt2,'-append','delimiter','')
        dlmwrite(file,txt3,'-append','delimiter','')
        dlmwrite(file,txt4,'-append','delimiter','')
        dlmwrite(file,txt5,'-append','delimiter','')
        dlmwrite(file,txt6,'-append','delimiter','')
        dlmwrite(file,txt7,'-append','delimiter','')
        
        %% output sce-file
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\scenario_output.sce',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('DIRECTORY("../outputlib/$1");');
        txt2 = sprintf('FILE("ElecProdSpec.out", "w")			= main.em.ep.opr.ElecProdSpec2;');
        txt3 = sprintf('FILE("GCap.out","w")                            = main.em.ep.gcp.GCap;');
        txt4 = sprintf('FILE("CO2Spec.out","w") 			= main.em.CO2Spec;');
        txt5 = sprintf('FILE("tpesEXT.out","w")                         = main.em.TPEStotalEXT;');
        txt6 = sprintf('FILE("TotalCostPerkWhNew.out", "w")		= main.em.ep.epc.TotalCostPerkWhNew;');
        txt7 = sprintf('FILE("TotalCostPerkWhAvg.out", "w")		= main.em.ep.epc.TotalCostPerkWhAvg;');
        txt8 = sprintf('FILE("StorLossISOTot.out", "w")			= main.em.ep.isp.StorLossISOTot;');
        txt9 = sprintf('FILE("StorLossFracISO.out", "w")		= main.em.ep.isp.StorLossFracISO;');
        txt10 = sprintf('FILE("ISODeple.out", "w")                       = main.em.ep.isod.ISODeple;');
        txt11 = sprintf('FILE("ISOAddOnNew.out", "w")			= main.em.ep.epc.ISOAddOnNew;');
        txt12 = sprintf('FILE("ISOAddOnAvg.out", "w")			= main.em.ep.epc.ISOAddOnAvg;');
        txt13 = sprintf('FILE("TotalCostperkWhAvgPVres.out","w")         = main.em.ep.epc.TotalCostperkWhAvgPVres;');
        txt14 = sprintf('FILE("BLF_Prod.out","w") 			= main.em.bio.BCap.BLFProdTot;');
        txt15 = sprintf('FILE("BSF_Prod.out","w") 			= main.em.bio.BCap.BSFProdTot;');
        txt16 = sprintf('FILE("tpes.out","w")                         	= main.em.TPEStotal;');
        
        dlmwrite(file,txt1,'')
        dlmwrite(file,txt2,'-append','delimiter','')
        dlmwrite(file,txt3,'-append','delimiter','')
        dlmwrite(file,txt4,'-append','delimiter','')
        dlmwrite(file,txt5,'-append','delimiter','')
        dlmwrite(file,txt6,'-append','delimiter','')
        dlmwrite(file,txt7,'-append','delimiter','')
        dlmwrite(file,txt8,'-append','delimiter','')
        dlmwrite(file,txt9,'-append','delimiter','')
        dlmwrite(file,txt10,'-append','delimiter','')
        dlmwrite(file,txt11,'-append','delimiter','')
        dlmwrite(file,txt12,'-append','delimiter','')
        dlmwrite(file,txt13,'-append','delimiter','')
        dlmwrite(file,txt14,'-append','delimiter','')
        dlmwrite(file,txt15,'-append','delimiter','')
        dlmwrite(file,txt16,'-append','delimiter','')
        
        %% settings bat-file
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\scenario_settings.bat',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('REM Scenarios settings');
        dlmwrite(file,txt1,'')
        
    end
end

%%

%% Time dependent output for TIMER without climate impacts, just historic values

clear tv jv CostCurveSmthOnshoret CostCurveSmthOffshoret LFCurveSmthOnshoret LFCurveSmthOffshoret MaxProdOnshoret MaxProdOffshoret
tv=[1971 2000 2050 2100];
% for historic only curves
jv=[1 1 1 1 1 1 1 1 ... % GFDL-ESM2M RCP26 hist hist hist; RCP60 hist hist hist
    6 6 6 6 6 6 6 6 ... % HadGEM-ES RCP26 hist hist hist; RCP60 hist hist hist
    11 11 11 11 11 11 11 11 ... %IPSL-CM5A-LR RCP26 hist hist hist; RCP60 hist hist hist
    16 16 16 16 16 16 16 16]; % MIROC5 RCP26 hist hist hist; RCP60 hist hist hist

j=0;
k=0;
l=0;
%Prep order curves
for i=1:numel(jv)
    j=j+1;
    if j>4; j=1; end
    if j==1; k=k+1; l=0; end
    %fprintf('t=%d j=%d jv=%d k=%d l=%d\n',i,j,jv(i),k,l)
    
    for c=1:101 %cost-curve steps
        for R=1:27 %region
            l=l+1;
            if c==1 && R==1
                CostCurveSmthOnshoret{k}(l)=tv(j);
                CostCurveSmthOffshoret{k}(l)=tv(j);
                LFCurveSmthOnshoret{k}(l)=tv(j);
                LFCurveSmthOffshoret{k}(l)=tv(j);
                l=l+1;
                CostCurveSmthOnshoret{k}(l)=CostCurveSmthOnshore{jv(i)}(c,R+1);
                CostCurveSmthOffshoret{k}(l)=CostCurveSmthOffshore{jv(i)}(c,R+1);
                LFCurveSmthOnshoret{k}(l)=LFCurveSmthOnshore{jv(i)}(c,R+1);
                LFCurveSmthOffshoret{k}(l)=LFCurveSmthOffshore{jv(i)}(c,R+1);
            else
                CostCurveSmthOnshoret{k}(l)=CostCurveSmthOnshore{jv(i)}(c,R+1);
                CostCurveSmthOffshoret{k}(l)=CostCurveSmthOffshore{jv(i)}(c,R+1);
                LFCurveSmthOnshoret{k}(l)=LFCurveSmthOnshore{jv(i)}(c,R+1);
                LFCurveSmthOffshoret{k}(l)=LFCurveSmthOffshore{jv(i)}(c,R+1);
            end
            
        end
        
    end
end

% prep order maxprod
j=0;
k=0;
l=0;
for i=1:numel(jv)
    j=j+1;
    if j>4; j=1; end
    if j==1; k=k+1; l=0; end
    %fprintf('t=%d j=%d jv=%d k=%d l=%d\n',i,j,jv(i),k,l)
    
    for R=1:27 %region
        l=l+1;
        if R==1
            MaxProdOnshoret{k}(l)=tv(j);
            MaxProdOffshoret{k}(l)=tv(j);
            l=l+1;
            MaxProdOnshoret{k}(l)=MaxProdOnshore{jv(i)}(R);
            MaxProdOffshoret{k}(l)=MaxProdOffshore{jv(i)}(R);
        else
            MaxProdOnshoret{k}(l)=MaxProdOnshore{jv(i)}(R);
            MaxProdOffshoret{k}(l)=MaxProdOffshore{jv(i)}(R);
        end
    end
end


% Writing
GCMIDt = {'GFDL-ESM2M' 'HADGEM2-ES' 'IPSL-CM5A-LR' 'MIROC5'};
RCPIDt = {'hist'};
c=0;
for i=1:4
    for j=1
        c=c+2;
        scenlib = sprintf('Y:\\Kennisbasis\\IMAGE\\model\\users\\david\\Pojects\\ISIMIP2E\\ISIMIP2E\\2_TIMER\\scenlib\\TIMER_2015\\ISIMIP2E\\Wind_%s_%s',GCMIDt{i},RCPIDt{j});
        matpath = fullfile(scenlib, sprintf('\\curves'));
        if ~isdir(matpath); mkdir(matpath); end
        pathname = fileparts(scenlib);
        
        % Onshore cost
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\CostCurveSmthWON.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.WONCostCurveSmth[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,CostCurveSmthOnshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % Offshore cost
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\CostCurveSmthWOFF.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.WOFFCostCurveSmth[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,CostCurveSmthOffshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % Onshore MaxProd
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\MaxProdWON.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.MaxProdOnshore[27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,MaxProdOnshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % Offshore MaxProd
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\MaxProdWOFF.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.MaxProdOffshore[27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,MaxProdOffshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % Onshore load
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\LoadCurveSmthWON.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.LFCurveSmthOnshore[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,LFCurveSmthOnshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        % Offshore load
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\curves\\LoadCurveSmthWOFF.dat',GCMIDt{i},RCPIDt{j}));
        txt=sprintf('real main.LFCurveSmthOffshore[101,27](t) = [');
        dlmwrite(file,txt,'');
        dlmwrite(file,LFCurveSmthOffshoret{c}','-append','delimiter','\t');
        dlmwrite(file,']','-append');
        dlmwrite(file,';','-append');
        
        %% input sce-file
        % Keep format. in the sce-file the format is correct
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\scenario_input.sce',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('DIRECTORY("../scenlib/$1/curves");');
        txt2 = sprintf('FILE("CostCurveSmthWON.dat","r")	=  main.em.ep.isod.CostCurveSmthWON;');
        txt3 = sprintf('FILE("CostCurveSmthWOFF.dat","r")	=  main.em.ep.isod.CostCurveSmthWOFF;');
        txt4 = sprintf('FILE("MaxProdWON.dat","r")              =  main.em.ep.isod.MaxProdWON;');
        txt5 = sprintf('FILE("MaxProdWOFF.dat","r")             =  main.em.ep.isod.MaxProdWOFF;');
        txt6 = sprintf('FILE("LoadCurveSmthWON.dat","r")	=  main.em.ep.isod.LoadCurveSmthWON;');
        txt7 = sprintf('FILE("LoadCurveSmthWOFF.dat","r")	=  main.em.ep.isod.LoadCurveSmthWOFF;');
        
        dlmwrite(file,txt1,'')
        dlmwrite(file,txt2,'-append','delimiter','')
        dlmwrite(file,txt3,'-append','delimiter','')
        dlmwrite(file,txt4,'-append','delimiter','')
        dlmwrite(file,txt5,'-append','delimiter','')
        dlmwrite(file,txt6,'-append','delimiter','')
        dlmwrite(file,txt7,'-append','delimiter','')
        
        %% output sce-file
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\scenario_output.sce',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('DIRECTORY("../outputlib/$1");');
        txt2 = sprintf('FILE("ElecProdSpec.out", "w")			= main.em.ep.opr.ElecProdSpec2;');
        txt3 = sprintf('FILE("GCap.out","w")                            = main.em.ep.gcp.GCap;');
        txt4 = sprintf('FILE("CO2Spec.out","w") 			= main.em.CO2Spec;');
        txt5 = sprintf('FILE("tpesEXT.out","w")                         = main.em.TPEStotalEXT;');
        txt6 = sprintf('FILE("TotalCostPerkWhNew.out", "w")		= main.em.ep.epc.TotalCostPerkWhNew;');
        txt7 = sprintf('FILE("TotalCostPerkWhAvg.out", "w")		= main.em.ep.epc.TotalCostPerkWhAvg;');
        txt8 = sprintf('FILE("StorLossISOTot.out", "w")			= main.em.ep.isp.StorLossISOTot;');
        txt9 = sprintf('FILE("StorLossFracISO.out", "w")		= main.em.ep.isp.StorLossFracISO;');
        txt10 = sprintf('FILE("ISODeple.out", "w")                       = main.em.ep.isod.ISODeple;');
        txt11 = sprintf('FILE("ISOAddOnNew.out", "w")			= main.em.ep.epc.ISOAddOnNew;');
        txt12 = sprintf('FILE("ISOAddOnAvg.out", "w")			= main.em.ep.epc.ISOAddOnAvg;');
        txt13 = sprintf('FILE("TotalCostperkWhAvgPVres.out","w")         = main.em.ep.epc.TotalCostperkWhAvgPVres;');
        txt14 = sprintf('FILE("BLF_Prod.out","w") 			= main.em.bio.BCap.BLFProdTot;');
        txt15 = sprintf('FILE("BSF_Prod.out","w") 			= main.em.bio.BCap.BSFProdTot;');
        txt16 = sprintf('FILE("tpes.out","w")                         	= main.em.TPEStotal;');
        
        dlmwrite(file,txt1,'')
        dlmwrite(file,txt2,'-append','delimiter','')
        dlmwrite(file,txt3,'-append','delimiter','')
        dlmwrite(file,txt4,'-append','delimiter','')
        dlmwrite(file,txt5,'-append','delimiter','')
        dlmwrite(file,txt6,'-append','delimiter','')
        dlmwrite(file,txt7,'-append','delimiter','')
        dlmwrite(file,txt8,'-append','delimiter','')
        dlmwrite(file,txt9,'-append','delimiter','')
        dlmwrite(file,txt10,'-append','delimiter','')
        dlmwrite(file,txt11,'-append','delimiter','')
        dlmwrite(file,txt12,'-append','delimiter','')
        dlmwrite(file,txt13,'-append','delimiter','')
        dlmwrite(file,txt14,'-append','delimiter','')
        dlmwrite(file,txt15,'-append','delimiter','')
        dlmwrite(file,txt16,'-append','delimiter','')
        
        %% settings bat-file
        file = fullfile(pathname, sprintf('\\Wind_%s_%s\\scenario_settings.bat',GCMIDt{i},RCPIDt{j}));
        txt1 = sprintf('REM Scenarios settings');
        dlmwrite(file,txt1,'')
        
    end
end

%% Visualization

% figure(1);clf;imagesc(CM1_data(:,:,1)'); axis image; colormap(jet)%colormap(parula)
%
% figure(2);clf;imagesc(CM1{13}(:,:)); axis image; colormap(jet)%colormap(parula)
%
% figure(3);clf;
% for m=1:12
%     subplot(3,4,m)
%     imagesc(CM1{m}); axis image; colormap(jet)
% end

