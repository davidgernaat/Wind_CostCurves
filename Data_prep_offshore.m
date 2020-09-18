%Data prep
clear all
root = 'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER\Wind_CC';

%% Distance2Land
disp('Calculate distance to land')
fname = sprintf('%s\\input\\imagemask_land.mat', root);
load(fname)
[nr nc] = size(imagemask);
% figure(1);clf;imagesc(imagemask);

%Create rasters with borders
win=50;
mindem=1;

%Create extra window for repmat
Dis_win = zeros((nr+2*win),(nc+2*win));
Land_win = zeros((nr+2*win),(nc+2*win),'int8');
[nrw,ncw]=size(Dis_win);

%Put land map in raster with window
counter=0;
for c=(win+1):(ncw-win);
    for r=(win+1):(nrw-win);
        counter=counter+1;
        Land_win(r,c)=imagemask(counter);
    end
end

% figure(2);clf; h = imagesc(Land_win); axis image; colormap(colordata); colorbar;

% Calculate distance
pold=0;
N=nr*nc;
i=0;

for r=(win+1):(nr+win)
    for c=(win+1):(nc+win)
        %NaNs for land
        if Land_win(r,c)==1; Dis_win(r,c)=NaN; continue; end
        
        % Build window around target cell
        firstrow = r-win;
        lastrow  = r+win;
        firstcol = c-win;
        lastcol  = c+win;
        
        Pwin = Land_win(firstrow:lastrow,firstcol:lastcol);
        
        % figure(3);clf;imagesc(Pwin);axis image;
        
        % Build repmat cell distance map
        rr = r-win;
        dx = repmat(abs(-win:win), [2*win+1, 1]);
        dx = dx;
        dy = dx';
        ds = sqrt(dx.^2+dy.^2);
        dsmax = max(ds);

        %Calc distance using repmat dis table
        near = Pwin==mindem;
        near_w = ds(near);
        dis = max(0.5,min(near_w))-0.5;
        if isempty(dis)==true; dis=dsmax(1)-0.5; end
        Dis_win(r,c)=dis;
        
    end
end

% check
% cmap=jet(256);
% cmap(1,:) = [1 1 1];
% figure(1);clf;imagesc(Dis_win);axis image;colormap(cmap);

% Cut correct window for output
OffDis = Dis_win(win+1:(nr+win),win+1:(nc+win)); % in cells
% figure(2);clf;imagesc(Dis_win);axis image;colormap(jet);

%% Bathymetry
disp('Bathymetry data')
fname = sprintf('%s\\input\\Offshore\\GRIDONE_1D.nc', root);
Z_prep = ncread(fname,'z');

nr=10801; %10801
nc=21601; %21601

z = zeros(nr,nc);
i=0;
for r=1:nr
    for c=1:nc
        i=i+1;
        z(r,c) = Z_prep(i);
    end
end

z = z(1:end-1,1:end-1);
% figure(1);clf;imagesc(z);axis image;

[nr nc] = size(z);

%Aggregate
scalingf = 30; %scaling factor 1min map to 0.5deg (30min)

rb = linspace(scalingf,nr,nr/scalingf);
cb = linspace(scalingf,nc,nc/scalingf);

% Create zeros map for output
B = zeros((size(z)/scalingf),'single');
[nrW, ncW] = size(B);

%
r1=1;
c1=1;
i=0;
for r=rb
    i=i+1;
    if i==1; r1=1; else r1=rb(i-1)+1; end
    r2=r;
    j=0;
    
    for c=cb
        j=j+1;
        if j==1; c1=1; else c1=cb(j-1)+1; end
        c2=c;
        clear BasinBlock MB
        BasinBlock = z(r1:r2,c1:c2);
        
        MB=mean(BasinBlock(:));
        
        B(i,j) = z(r1,c1);
        
    end
end

% figure(1);clf;imagesc(B);axis image; colormap(jet); colorbar
% 
% rr=75; %Dutch coastal cell
% cc=369;
% B(rr,cc)

%
% cmap=jet(256);
% cmap(1,:) = [1 1 1];
% Bs=B;
% Bs(find(imagemask))=NaN;
% Bl=B;
% Bl(find(~imagemask))=NaN;
% figure(1);clf;
% ax1=subplot(1,2,1);imagesc(Bs);axis image; colormap(cmap); colorbar; title('B sea');
% ax2=subplot(1,2,2);imagesc(Bl);axis image; colormap(cmap); colorbar; title('B land');
% linkaxes([ax1 ax2])

%% Ice data
disp('Ice data')
fname = sprintf('%s\\input\\Offshore\\ICE_COVER.dat', root);
[Ice_data] = read_mym(fname);

Ice = zeros(360,720);
[nr nc] = size(Ice);

i=0;
for r=1:nr
    for c=1:nc
        i=i+1;
        Ice(r,c)=~Ice_data(i);
        if imagemask(r,c)==1; Ice(r,c)=1; end
    end
end

% figure(1);clf;imagesc(Ice);axis image

%% Shipping data
disp('Shipping data')
fname = sprintf('%s\\input\\Offshore\\Shipping_MyM_file.dat', root);
[Shipping_data] = read_mym(fname);

Ship = ones(360,720);
[nr nc] = size(Ship);

i=0;
for r=1:nr
    for c=1:nc
        i=i+1;
        if Shipping_data(i)==1
            Ship(r,c)=Shipping_data(i) * 0.8;
        end
    end
end

% figure(1);clf;imagesc(Ship);axis image

%% Marine protected data
disp('Marine protected areas')
fname = sprintf('%s\\input\\Offshore\\MPA_MyM_file.dat', root);
[MP_data] = read_mym(fname);

MP = ones(360,720);
[nr nc] = size(MP);

i=0;
for r=1:nr
    for c=1:nc
        i=i+1;
        if MP_data(i)==1
            MP(r,c)=0;
        end
    end
end

% figure(1);clf;imagesc(MP);axis image

%% Wind speeds NASA (m/s)
disp('Winds speeds NASA')
fname = sprintf('%s\\input\\offshore\\10yr_wspd50m.dat', root);
wsp_data = dlmread(fname,' ',7,0);

latv = linspace(-90,89,180/1);
lonv = linspace(-180,179,360/1);
for m=1:13
    wsp{m} = zeros(numel(latv),numel(lonv));
    [nr, nc] = size(wsp{m});
    for r=1:nr
        for c=1:nc
            
            [rlat, ~] = ind2sub(size(wsp_data),find(wsp_data(:,1)==latv(r)));
            [rlon, ~] = ind2sub(size(wsp_data),find(wsp_data(:,2)==lonv(c)));
            
            if sum(ismember(rlon,rlat)) > 0
                
                wsp{m}(r,c) = wsp_data(rlon(find(ismember(rlon,rlat))),m+2);
                
            end
            
        end
    end
    wsp{m} = flipud(wsp{m});
end

% figure(1);clf;
% for m=1:12
%     subplot(3,4,m); imagesc(wsp{m});colormap(jet)
% end

%% Wind speed difference 10 or 50m NASA (%)
disp('Wind speeds difference 10 or 50m')
fname = sprintf('%s\\input\\offshore\\10yr_pct10m_wnd.dat', root);
pct10m_data = dlmread(fname,' ',9,0);

latv = linspace(-90,89,180/1);
lonv = linspace(-180,179,360/1);
for m=1:13
    pct10m{m} = zeros(numel(latv),numel(lonv));
    [nr, nc] = size(pct10m{m});
    for r=1:nr
        for c=1:nc
            
            [rlat, ~] = ind2sub(size(pct10m_data),find(pct10m_data(:,1)==latv(r)));
            [rlon, ~] = ind2sub(size(pct10m_data),find(pct10m_data(:,2)==lonv(c)));
            
            if sum(ismember(rlon,rlat)) > 0
                
                pct10m{m}(r,c) = pct10m_data(rlon(find(ismember(rlon,rlat))),m+2);
                
            end
            
        end
    end
    pct10m{m} = flipud(pct10m{m});
end

% figure(2);clf; imagesc(pct10m{13});colormap(jet);colorbar

%% EEZ zones data
disp('EEZ zones')
fname = sprintf('%s\\input\\Offshore\\TIMER_REGION_REGIO28.dat', root);
EEZ_data = dlmread(fname,'',1,0);

EEZ = zeros(360,720);
[nr nc] = size(EEZ);

i=0;
for r=1:nr
    for c=1:nc
        i=i+1;
        if imagemask(r,c)==1; EEZ(r,c)=0; continue; end
        EEZ(r,c)=EEZ_data(i);
    end
end

% Bathymetry data correction, to removed cells that are inlands
[nr nc]=size(EEZ);
for r=1:nr
    for c=1:nc
        if EEZ(r,c)>0 && B(r,c)>0; B(r,c)=0; end
    end
end

% figure(3);clf;imagesc(B);axis image

%% EEZ zones countries
fname = sprintf('%s//input//Offshore//EEZ//EEZ.tif',root);

CEEZ_data = imread(fname);

% figure(2);clf;imagesc(CEEZ_data);colormap(prism)

% Georeference of shape file EEZ_maritimeboundaries_2009_v5_VLIZ
% 87,02
% -78,57
CEEZ = zeros(360,720);
CEEZ(7:end-23,:) = CEEZ_data;
CEEZ(CEEZ(:)==255) = 0;

[nr,nc] = size(CEEZ);
i=0;
for r=1:nr
    for c=1:nc
        i=i+1;
        if imagemask(r,c)==1; CEEZ(r,c)=0; end
    end
end

% Now match these ObjectID numbers to ISOcodes that are used in the rest of the code
% Read data from the ShapeFile attriute table, already matched with iso code in excel
fname = sprintf('%s\\input\\Offshore\\EEZ\\EEZ_codes.csv', root);
fileID = fopen(fname);
CC = textscan(fileID,'%s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',',','HeaderLines',1);
fclose(fileID);
for i=1:numel(CC{1})
    ObjID(i,1) = str2num(CC{1}{i})'; %ObjectID number
    ObjID(i,2) = str2num(CC{13}{i})'; %ObjectID number
end

i=0;
for r=1:nr
    for c=1:nc
        i=i+1;
        if CEEZ(r,c)==0; continue; end;
        CEEZ(r,c) = ObjID(find(ObjID(:,1)==CEEZ(r,c)),2);
    end
end

% figure(1);clf;
% ax1=subplot(2,1,1);
% imagesc(EEZ);colormap(jet);axis image
% ax2=subplot(2,1,2);
% imagesc(CEEZ);colormap(jet);axis image
% linkaxes([ax1 ax2])

%% save
disp('Save')
matfile = fullfile(root, sprintf('input\\input_data_offshore.mat'));
save(matfile,'OffDis','B','Ice','Ship','MP','wsp','pct10m','EEZ','CEEZ');

matfile = fullfile(root, sprintf('input\\input_data_offshore_ISIMIP.mat'));
save(matfile,'OffDis','B','Ice','Ship','MP','pct10m','EEZ','CEEZ');
