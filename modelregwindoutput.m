function modelregoutput(root,C2R_fname,scenlib,CMc,GCMID2,RCPID2,TIMEID2)

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

%% Average wind speed onshore

windspeedonshore = zeros(1,numel(unique_regs));
windspeedoffshore = zeros(1,numel(unique_regs));

for j=1:numel(unique_regs)
    windspeedonshore(:,j) = mean(CMc{13}(CRinda{j})); %m/s
    
    if isempty(CRindao{j})==1; windspeedoffshore(:,j)=0; continue; end;
    windspeedoffshore(:,j) = mean(CMc{13}(CRindao{j})); %m/s
end
% figure(3);clf;imagesc(CMc{13});axis image; colorbar

% Txt for output 
c=0;
for j=1:numel(unique_regs)
    if j==1
       c=c+1;
       txt1{c}=sprintf('unit: mean m/s over cells in region | Column ',j, unique_regs{j});
        
    end
    c=c+1;
    txt1{c}=sprintf('%d=%s;',j, unique_regs{j});
end
txt = horzcat(txt1{:});

% onshore
file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\Mean_Windspeed_%s.dat',GCMID2,RCPID2,TIMEID2,'onshore'));
dlmwrite(file,txt,'');
dlmwrite(file,round(windspeedonshore,3),'-append','delimiter',';');

% onshore
file = fullfile(pathname, sprintf('\\Wind_%s_%s_%s\\Mean_Windspeed_%s.dat',GCMID2,RCPID2,TIMEID2,'offshore'));
dlmwrite(file,txt,'');
dlmwrite(file,round(windspeedoffshore,3),'-append','delimiter',';');

