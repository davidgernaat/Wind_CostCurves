%% Run ISIMIPs
clear all

root = 'Y:\ontwapps\Timer\Users\David\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER\Wind_CC';
root_data = 'Y:\ontwapps\Timer\Users\David\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER\Wind_CC\input\ISIMIP';

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
        fprintf('%s\n',fnames(j).name(13:70))
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
% i=3;

for i=1:numel(names)
    fprintf('Reading file #%d of %d\n',i, numel(names))
    clear CM 
    
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
                CMc{i}{m}(r,c) = CM{m}(r,c) * 1;
            end
        end
    end
end

figure(1);clf;imagesc(CMc{1}{13});colorbar

%Aswan dam cell Egypt
r=130;
c=422;

for ne=1:numel(CMc)
    fprintf('%d %s %s %s %0.2f unit\n',ne,GCMID{ne},RCPID{ne},TIMEID{ne},CMc{ne}{13}(r,c))
end

