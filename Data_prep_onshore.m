%Data prep
clear all
root = 'Y:\Kennisbasis\IMAGE\model\users\david\Pojects\ISIMIP2E\ISIMIP2E\2_TIMER\Wind_CC';

%% Image land mask
disp('IMAGE land mask')
fname = sprintf('%s\\input\\imagemask_land.mat', root);
load(fname)

[nr nc] = size(imagemask);

%% Read annual wind speed data CRU database - this section is blocked because it is calculated one section below
% disp('Annual wind speed data')
% fname = sprintf('%s\\input\\main.vavg10.dat', root);
% [V_data, t] = read_mym(fname);
% 
% V = zeros(360,720);
% 
% i=0;
% for r=1:nr
%     for c=1:nc
%         if imagemask(r,c)==0; continue; end;
%         i=i+1;
%         V(r,c)=V_data(i);
%     end
% end
% figure(1);clf;imagesc(V);axis image; colormap(jet); colorbar

%% Read monthly wind speed data CRU database (m/s)
disp('Monthly wind speed data')
fname = sprintf('%s\\input\\main.GWND6190.dat', root);
[Vm_data, t] = read_mym(fname);

for m=1:12
    V{m} = zeros(360,720);
    i=0;
    for r=1:nr
        for c=1:nc
            if imagemask(r,c)==0; continue; end;
            i=i+1;
            V{m}(r,c)=Vm_data(1,m,i);
        end
    end
end

% Yearly average
V{13} = zeros(360,720);
for r=1:nr
    for c=1:nc
        for m=1:12
            dt(m) = V{m}(r,c);
        end
        dt = mean(dt);
        V{13}(r,c) = dt;
    end
end

% figure(1);clf;
% for m=1:12
%     subplot(4,3,m)
%     imagesc(V{m});axis image; colormap(jet)
% end
% 
% figure(1);clf;imagesc(V{13});axis image; colormap(jet); colorbar

%%
disp('Altitude')
fname = sprintf('%s\\input\\main.ALT.dat', root);
[Alt_data t] = read_mym(fname);

Alt = zeros(360,720);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        Alt(r,c)=Alt_data(i);
    end
end

% figure(1);clf;imagesc(Alt);axis image; colormap(jet)

%% IMAGE regions
disp('IMAGE Regions')
fname = sprintf('%s\\input\\region27.dat', root);
[Region_data] = read_mym(fname);

IRegion = zeros(360,720);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        IRegion(r,c)=Region_data(i);
    end
end

% figure(1);clf;imagesc(IRegion);axis image; colormap(jet)

%% Countries
disp('Countries')
fname = sprintf('%s\\input\\gISO2.asc', root);
GISO = zeros(360,720);
GISO(2:end,2:end) = dlmread(fname,',',7,1);
GISO(GISO==-9999)=NaN;
% imagesc(GISO); colorbar

%% Area regions
disp('Area per cell')
Area = zeros(360,720); %km2
latv = linspace(-89.5,89.5,180/0.5);

for r=1:nr
    for c=1:nc
        Area(r,c)=55*(55*cos(deg2rad(latv(r))));
    end
end
% r=171;
% c=212;
% r=13; %first imagemask row number compatible with MyM
% c=284; %first imagemask column number compatible with MyM
% figure(1);clf;imagesc(Area);axis image; colormap(jet); colorbar;hold on
% plot(c,r,'.k','markersize',20)
% hold off

fname = sprintf('%s\\input\\area.dat', root);
[Area_data] = read_mym(fname);

AreaHoogwijk = zeros(360,720); %Landcover 2010

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        AreaHoogwijk(r,c)=Area_data(1,i); %Seems to be incorrect in the high latitude cells, use the self-calculated one
    end
end

% figure(2);clf;imagesc(AreaHoogwijk);axis image; colormap(jet); colorbar; hold on
% plot(c,r,'.k','markersize',20)
% hold off

%% GLCT
disp('GLCT')
fname = sprintf('%s\\input\\GLCT_SSP2.asc', root);
[GLCT_data t] = read_mym(fname);

% Put in map
GLCT = zeros(360,720); %Landcover 2010

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        GLCT(r,c)=GLCT_data(9,i);
    end
end

% figure(1);clf;imagesc(GLCT);axis image; colormap(jet); colorbar

%% GLCT data from winddat --> wind_27
disp('GLCT wind_27 winddat')
fname = sprintf('%s\\input\\glct_winddat.out', root);
[glctwindat_data t] = read_mym(fname);

glctwindat = zeros(360,720);
[nr nc] = size(glctwindat);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        glctwindat(r,c)=glctwindat_data(5,i);
    end
end

% figure(1);clf;imagesc(glctwindat);axis image

%%
disp('Abandoned agricultural land wind_27 winddat')
fname = sprintf('%s\\input\\AbonAgLand.out', root);
[abonAgLand_data t] = read_mym(fname);

abonAgLand = zeros(360,720);
[nr nc] = size(abonAgLand);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        abonAgLand(r,c)=abonAgLand_data(5,i);
    end
end

% figure(1);clf;imagesc(abonAgLand);axis image

%% Bioreservce
disp('Bio reserves')
fname = sprintf('%s\\input\\BIORESERVE.dat', root);
[BIOres_data t] = read_mym(fname);

% Put in map
BIOres = zeros(360,720); %Bioreserves [type 0, 1 or 2], value of 0 means available

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        BIOres(r,c)=BIOres_data(i);
    end
end

% figure(1);clf;imagesc(BIOres);axis image; colormap(jet); colorbar

%% Buildup
disp('Build Up')
fname = sprintf('%s\\input\\BUILDUP.dat', root);
[Buildup_data t] = read_mym(fname);

% Put in map
Buildup = zeros(360,720); %Urban areas (Fraction per cell)

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        Buildup(r,c)=Buildup_data(i);
    end
end

% figure(1);clf;imagesc(Buildup);axis image; colormap(jet); colorbar

%% Pop
disp('Population')
fname = sprintf('%s\\input\\totpop.dat', root);
[totpop_data t] = read_mym(fname);
fname = sprintf('%s\\input\\urbpop.dat', root);
[urbpop_data t] = read_mym(fname);
fname = sprintf('%s\\input\\rurpop.dat', root);
[rurpop_data t] = read_mym(fname);

% Put in map
totpop = zeros(360,720);
urbpop = zeros(360,720);
rurpop = zeros(360,720);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        totpop(r,c)=totpop_data(i);
        urbpop(r,c)=urbpop_data(i);
        rurpop(r,c)=rurpop_data(i);
    end
end

% figure(1);clf;imagesc(log(totpop));axis image; colormap(jet); colorbar
% figure(2);clf;imagesc(log(urbpop));axis image; colormap(jet); colorbar
% figure(3);clf;imagesc(log(rurpop));axis image; colormap(jet); colorbar

%% Slope
disp('Slope')
fname = sprintf('%s\\input\\Slope_lessthan_3pcnt.dat', root);
[Slope_data t] = read_mym(fname);

% Put in map
Slope = zeros(360,720); %!Fraction of cell area with slope <3%

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        Slope(r,c)=Slope_data(i);
    end
end

% figure(1);clf;imagesc(Slope);axis image; colormap(jet); colorbar

%% DistanceToLoad
disp('Distance to load')

fname = sprintf('%s\\input\\transmission2010.dat', root);
[DistanceToLoad_data ~] = read_mym(fname);

% Put in map
DistanceToLoad = zeros(360,720);

i=0;
for r=1:nr
    for c=1:nc
        if imagemask(r,c)==0; continue; end;
        i=i+1;
        if DistanceToLoad_data(1,i)==0; DistanceToLoad(r,c)=1; continue; end;
        DistanceToLoad(r,c)=DistanceToLoad_data(1,i);
    end
end

% DistanceToLoad(DistanceToLoad(:)>1000)=1000;

% figure(1);clf;imagesc(log(totpop));axis image; colormap(jet); colorbar

%% IMAGE regions names

IMAGER{1}='Canada';
IMAGER{2}='USA';
IMAGER{3}='Mexico';
IMAGER{4}='Rest Central America';
IMAGER{5}='Brazil';
IMAGER{6}='Rest South America';
IMAGER{7}='Northern Africa';
IMAGER{8}='Western Africa';
IMAGER{9}='Eastern Africa';
IMAGER{10}='South Africa';
IMAGER{11}='Western Europe';
IMAGER{12}='Central Europe';
IMAGER{13}='Turkey';
IMAGER{14}='Ukraine +';
IMAGER{15}='Asia-Stan';
IMAGER{16}='Russia +';
IMAGER{17}='Middle East';
IMAGER{18}='India +';
IMAGER{19}='Korea';
IMAGER{20}='China +';
IMAGER{21}='Southeastern Asia';
IMAGER{22}='Indonesia +';
IMAGER{23}='Japan';
IMAGER{24}='Oceania';
IMAGER{25}='Rest S.Asia';
IMAGER{26}='Rest S.Africa';
IMAGER{27}='World';

%% Reading GDPpc ppp
fname = sprintf('%s\\input\\GDPpc2010IsoCode.csv', root);
fileID = fopen(fname);
C = textscan(fileID,'%s %s %s %s %s %s %s','Delimiter',';','HeaderLines',1);
fclose(fileID);
for i=1:numel(C{3})
    ISOGDP(i,1) = str2num(C{3}{i})'; %ISO number
    ISOGDP(i,2) = str2num(C{4}{i})'; %GDP pc
    ISOGDP(i,3) = str2num(C{5}{i})'; %GDP MER
    
    CountryNames{i} = C{1}{i};       %Country names
end

%% save
disp('Save')

matfile = fullfile(root, sprintf('input\\input_data_onshore.mat'));
save(matfile,'IRegion','Area','GLCT','BIOres','Buildup','imagemask','totpop','urbpop','rurpop','Slope','DistanceToLoad','IMAGER', 'V', 'Alt','glctwindat','abonAgLand','GISO','AreaHoogwijk','ISOGDP','CountryNames');

matfile = fullfile(root, sprintf('input\\input_data_onshore_ISIMIP.mat'));
save(matfile,'IRegion','Area','GLCT','BIOres','Buildup','imagemask','totpop','urbpop','rurpop','Slope','DistanceToLoad','IMAGER', 'Alt','glctwindat','abonAgLand','GISO','AreaHoogwijk','ISOGDP','CountryNames');
