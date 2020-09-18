function modelregoutput(root,C2R_fname,scenlib,AnnTechPotCellCSP2,CSP_COE_PL2,GCMID2,RCPID2,TIMEID2)

% i=6;
% root;
% C2R_fname = sprintf('%s\\input\\ISIMIP\\Modelregionallocation\\Country_to_POLESregion_IsoCode.csv',root);
% scenlib =sprintf('%s\\scenlib\\TIMER_2015\\ISIMIP2E\\POLES\\Wind_%s_%s_%s',root2,GCMID{i},RCPID{i},TIMEID{i});
% AnnTechPotCellCSP2 = TechPotCell{i}{13};
% CSP_COE_PL2 = COE{i}{13};
% GCMID2 = GCMID{i};
% RCPID2 = RCPID{i};
% TIMEID2 = TIMEID{i};

fname = sprintf('%s\\input\\input_data_onshore.mat', root);
load(fname)

fname = sprintf('%s\\input\\input_data_offshore.mat', root);
load(fname)

%% Model region names

% C2R_fname=sprintf('%s\\data\\ISIMIP\\Modelregionallocation\\Country_to_REMINDregion_IsoCode.csv',root);

fid   = fopen(C2R_fname);
country_reg = textscan(fid, '%s %s %s %s','Delimiter',',','HeaderLines',1);
fclose(fid);

for j=1:numel(country_reg{3})
    ISON(j,1) = str2num(country_reg{3}{j})'; %ISO number countries
end

unique_regs = unique(country_reg{4}); %model region names

%% find model region index numbers

% Model regs indices onshore
for j=1:numel(unique_regs)
    modelregs_ISONS{j} = ISON(find(strcmp(country_reg{4}, unique_regs{j})),1); %Which country ISO numbers belong to the model region

    for k=1:numel(modelregs_ISONS{j})
        CRind{j}{k} = find(GISO(:)==modelregs_ISONS{j}(k));
    end
    CRinda{j} = unique(vertcat(CRind{j}{:}));
    
end


% Model regs indices offshore
for j=1:numel(unique_regs)
    
    for k=1:numel(modelregs_ISONS{j})
        CRindo{j}{k} = find(CEEZ(:)==modelregs_ISONS{j}(k));
    end
    CRindao{j} = unique(vertcat(CRindo{j}{:}));
    
end

%% check
% dt=CEEZ;
% dt(dt>1)=1;
% dt(dt==0)=1;
% dt(isnan(dt))=-1;
% dt(CRindao{1})=2;
% figure(2);clf;imagesc(dt);axis image; colorbar

%%
% scenlib = sprintf('Y:\\Kennisbasis\\IMAGE\\model\\Users\\David\\Pojects\\ISIMIP2E\\ISIMIP2E\\2_TIMER\\scenlib\\TIMER_2015\\ISIMIP2E\\REMIND\\Solar_%s_%s_%s',GCMID{i},RCPID{i},TIMEID{i});

matpath = fullfile(scenlib, sprintf(''));
if ~isdir(matpath)
    mkdir(matpath);
end

pathname = fileparts(scenlib);

%% Maxprods onshore

MaxProdCSP = zeros(1,numel(unique_regs)+1);

for j=1:numel(unique_regs)
    MaxProdCSP(:,j) = sum(AnnTechPotCellCSP2(CRinda{j})); %kWh 
    MaxProdCSPo(:,j) = sum(AnnTechPotCellCSP2(CRindao{j})); %kWh 
end
MaxProdCSP(:,numel(unique_regs)+1) =sum(MaxProdCSP(:,1:numel(unique_regs)));
MaxProdCSPo(:,numel(unique_regs)+1) =sum(MaxProdCSPo(:,1:numel(unique_regs)));

% figure(3);clf;imagesc(AnnTechPotCellPV);axis image; colorbar

% Txt for output maxprod
c=0;
for j=1:numel(unique_regs)
    if j==1
       c=c+1;
       txt1{c}=sprintf('unit: kWh/y | Column ',j, unique_regs{j});
        
    end
    c=c+1;
    txt1{c}=sprintf('%d=%s;',j, unique_regs{j});
end
txt = horzcat(txt1{:});

% onshore
file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\MaxProd_%s.dat',GCMID2,RCPID2,TIMEID2,'onshore'));
dlmwrite(file,txt,'');
dlmwrite(file,round(MaxProdCSP,3),'-append','delimiter',';');

% onshore
file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\MaxProd_%s.dat',GCMID2,RCPID2,TIMEID2,'offshore'));
dlmwrite(file,txt,'');
dlmwrite(file,round(MaxProdCSPo,3),'-append','delimiter',';');

%% CostCurves prep 
nrs=numel(unique_regs)+1;

% 
clear SI
for j=1:numel(unique_regs)
    CSP_COE_IR{j} = round(CSP_COE_PL2(CRinda{j}),3);
    CSP_Pot_IR{j} = AnnTechPotCellCSP2(CRinda{j});
    
    CSP_COE_IRo{j} = round(CSP_COE_PL2(CRindao{j}),3);
    CSP_Pot_IRo{j} = AnnTechPotCellCSP2(CRindao{j});
end

CSP_COE_IR{numel(unique_regs)+1} = CSP_COE_PL2(find(isnan(GISO)));
CSP_Pot_IR{numel(unique_regs)+1} = AnnTechPotCellCSP2(find(isnan(GISO)));

CSP_COE_IRo{numel(unique_regs)+1} = CSP_COE_PL2(CEEZ>0);
CSP_Pot_IRo{numel(unique_regs)+1} = AnnTechPotCellCSP2(CEEZ>0);

for j=1:numel(unique_regs)
    [CSP_COE_IR_S{j}, SI{j}] = sort(CSP_COE_IR{j});
    CSP_Pot_IR_S{j} = CSP_Pot_IR{j}(SI{j});
    
    [CSP_COE_IR_So{j}, SIo{j}] = sort(CSP_COE_IRo{j});
    CSP_Pot_IR_So{j} = CSP_Pot_IRo{j}(SIo{j});
end

[CSP_COE_IR_S{nrs}, SI{nrs}] = sort(CSP_COE_IR{nrs});
CSP_Pot_IR_S{nrs} = CSP_Pot_IR{nrs}(SI{nrs});

[CSP_COE_IR_So{nrs}, SIo{nrs}] = sort(CSP_COE_IRo{nrs});
CSP_Pot_IR_So{nrs} = CSP_Pot_IRo{nrs}(SIo{nrs});

for j=1:nrs
    CSP_COE_IR_S{j} = CSP_COE_IR_S{j}(~isnan(CSP_COE_IR_S{j}));
    CSP_Pot_IR_S{j} = CSP_Pot_IR_S{j}(~isnan(CSP_COE_IR_S{j}));
    
    CSP_COE_IR_S{j} = CSP_COE_IR_S{j}(find(~isinf(CSP_COE_IR_S{j})));
    CSP_Pot_IR_S{j} = CSP_Pot_IR_S{j}(find(~isinf(CSP_COE_IR_S{j})));
    
    CSP_COE_IR_So{j} = CSP_COE_IR_So{j}(~isnan(CSP_COE_IR_So{j}));
    CSP_Pot_IR_So{j} = CSP_Pot_IR_So{j}(~isnan(CSP_COE_IR_So{j}));
    
    CSP_COE_IR_So{j} = CSP_COE_IR_So{j}(find(~isinf(CSP_COE_IR_So{j})));
    CSP_Pot_IR_So{j} = CSP_Pot_IR_So{j}(find(~isinf(CSP_COE_IR_So{j})));
end

%% Transform into TIMER curves

for j=1:nrs
    
    %onshore
    CSP_Pot_cumsum{j} = cumsum(CSP_Pot_IR_S{j}); %kWh
    
    for i=1:numel(CSP_Pot_cumsum{j})
        CSP_Pot_Indexed{j}(i) = CSP_Pot_cumsum{j}(i)/CSP_Pot_cumsum{j}(end); 
    end
    
    %offshore
    CSP_Pot_cumsumo{j} = cumsum(CSP_Pot_IR_So{j}); %kWh
    
    for i=1:numel(CSP_Pot_cumsumo{j})
        CSP_Pot_Indexedo{j}(i) = CSP_Pot_cumsumo{j}(i)/CSP_Pot_cumsumo{j}(end); 
    end
    
end

%% Find min idx
val = linspace(0,1,100);

for j=1:nrs
    for i=1:numel(val)
        
        %onshore
        clear tmp
        if isempty(CSP_Pot_Indexed{j})==1; continue; end
        tmp = abs((CSP_Pot_Indexed{j}./val(i))-1);
        [~, idxminCSP{j}(i)] = min(tmp); %index of closest value
        
        %offshore
        clear tmp
        if isempty(CSP_Pot_Indexedo{j})==1; continue; end
        tmp = abs((CSP_Pot_Indexedo{j}./val(i))-1);
        [~, idxminCSPo{j}(i)] = min(tmp); %index of closest value
        
    end
end

%% Put in TIMER matrix

a = [0	0.01	0.02	0.03	0.04	0.05	0.06	0.07	0.08	0.09	0.1	0.11	0.12	0.13	0.14	0.15	0.16	0.17	0.18	0.19	0.2	0.21	0.22	0.23	0.24	0.25	0.26	0.27	0.28	0.29	0.3	0.31	0.32	0.33	0.34	0.35	0.36	0.37	0.38	0.39	0.4	0.41	0.42	0.43	0.44	0.45	0.46	0.47	0.48	0.49	0.5	0.51	0.52	0.53	0.54	0.55	0.56	0.57	0.58	0.59	0.6	0.61	0.62	0.63	0.64	0.65	0.66	0.67	0.68	0.69	0.7	0.71	0.72	0.73	0.74	0.75	0.76	0.77	0.78	0.79	0.8	0.81	0.82	0.83	0.84	0.85	0.86	0.87	0.88	0.89	0.9	0.91	0.92	0.93	0.94	0.95	0.96	0.97	0.98	0.99	1];

CostCurveSmthCSP = zeros(100,nrs+1);
CostCurveSmthCSPo = zeros(100,nrs+1);

for j=1:numel(a)-1
    CostCurveSmthCSP(j,1) = a(j);
    CostCurveSmthCSPo(j,1) = a(j);
end

for j=1:nrs
    COEISort_Indexed_CSP{j} = CSP_COE_IR_S{j}(idxminCSP{j});
    COEISort_Indexed_CSPo{j} = CSP_COE_IR_So{j}(idxminCSPo{j});
end

for j=1:nrs
    % Cost curves
    % onshore
    if isempty(COEISort_Indexed_CSP{j})==1; continue; end
    CostCurveSmthCSP(:,j+1) = COEISort_Indexed_CSP{j};
    
    % offshore
    if isempty(COEISort_Indexed_CSPo{j})==1; continue; end
    CostCurveSmthCSPo(:,j+1) = COEISort_Indexed_CSPo{j};
    
end

for j=1:nrs
    % Cost curves
    % onshore
    CostCurveSmthCSP(101,1) = 1;
    
    if isempty(COEISort_Indexed_CSP{j})==1
        CostCurveSmthCSP(101,j+1) = 0;
    else
        CostCurveSmthCSP(101,j+1) = COEISort_Indexed_CSP{j}(end);
    end
    
    % offshore
    CostCurveSmthCSPo(101,1) = 1;
    
    if isempty(COEISort_Indexed_CSPo{j})==1
        CostCurveSmthCSPo(101,j+1) = 0;
    else
        CostCurveSmthCSPo(101,j+1) = COEISort_Indexed_CSPo{j}(end);
    end
end

%% Txt for output costcurve
c=0;
for j=1:numel(unique_regs)
    if j==1
       c=c+1;
       txt1{c}=sprintf('unit: $/kWh | Column 1=fraction of maximum potential (maxprod);',j, unique_regs{j});
        
    end
    c=c+1;
    txt1{c}=sprintf('%d=%s;',j+1, unique_regs{j});
end
txt = horzcat(txt1{:});

% onshore
file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\CostCurve_%s.dat',GCMID2,RCPID2,TIMEID2,'onshore'));
dlmwrite(file,txt,'');
dlmwrite(file,round(CostCurveSmthCSP,3),'-append','delimiter',';');

% offshore
file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\CostCurve_%s.dat',GCMID2,RCPID2,TIMEID2,'offshore'));
dlmwrite(file,txt,'');
dlmwrite(file,round(CostCurveSmthCSPo,3),'-append','delimiter',';');
