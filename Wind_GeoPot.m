
%% Suitability factors

% if abandoned agricultural land than suitability factor 0.9

% !1	agricultural land Hoogwijk 
% !2	extensive grassland
% !3	carbon plantations
% !4	regrowth forest (Abandoning)
% !5	regrowth forest (Timber)
% !6	biofuel
% !7	ice
% !8	tundra 
% !9	wooded tundra
% !10	boreal forest
% !11	cool conifer
% !12	temp. mixed forest
% !13	temp.deciduous forest
% !14	warm mixed
% !15	grassland/steppe 
% !16	hot desert 
% !17	shrubland 
% !18	savanna 
% !19	tropical woodland
% !20	tropical forest

if glctfile==2 %new glct file
    SuitabilityFactor = ones(size(GLCT)); %ones for offshore
    for r=1:nr
        for c=1:nc
            if GLCT(r,c) == 1; SuitabilityFactor(r,c) = 0.7;
            elseif GLCT(r,c) == 2; SuitabilityFactor(r,c) = 0.8;
            elseif GLCT(r,c) == 3; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 4; SuitabilityFactor(r,c) = 0.1;
            elseif GLCT(r,c) == 5; SuitabilityFactor(r,c) = 0.1;
            elseif GLCT(r,c) == 6; SuitabilityFactor(r,c) = 0.7;
            elseif GLCT(r,c) == 7; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 8; SuitabilityFactor(r,c) = 0.8;
            elseif GLCT(r,c) == 9; SuitabilityFactor(r,c) = 0.5;
            elseif GLCT(r,c) == 10; SuitabilityFactor(r,c) = 0.1;
            elseif GLCT(r,c) == 11; SuitabilityFactor(r,c) = 0.1;
            elseif GLCT(r,c) == 12; SuitabilityFactor(r,c) = 0.1;
            elseif GLCT(r,c) == 13; SuitabilityFactor(r,c) = 0.1;
            elseif GLCT(r,c) == 14; SuitabilityFactor(r,c) = 0.1;
            elseif GLCT(r,c) == 15; SuitabilityFactor(r,c) = 0.8;
            elseif GLCT(r,c) == 16; SuitabilityFactor(r,c) = 1;
            elseif GLCT(r,c) == 17; SuitabilityFactor(r,c) = 0.5;
            elseif GLCT(r,c) == 18; SuitabilityFactor(r,c) = 0.9;
            elseif GLCT(r,c) == 19; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 20; SuitabilityFactor(r,c) = 0;
            end
            
            if abonAgLand(r,c) ==1
                GLCT(r,c) = 0.8;
            end
        end
    end
else %hoogwijk glct file
    
    SuitabilityFactor = ones(size(GLCT)); %ones for offshore
    for r=1:nr
        for c=1:nc
            if GLCT(r,c) == 1; SuitabilityFactor(r,c) = 0.5;
            elseif GLCT(r,c) == 2; SuitabilityFactor(r,c) = 0.25;
            elseif GLCT(r,c) == 3; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 4; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 5; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 6; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 7; SuitabilityFactor(r,c) = 0.25;
            elseif GLCT(r,c) == 8; SuitabilityFactor(r,c) = 0.25;
            elseif GLCT(r,c) == 9; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 10; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 11; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 12; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 13; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 14; SuitabilityFactor(r,c) = 0.25;
            elseif GLCT(r,c) == 15; SuitabilityFactor(r,c) = 0.9;
            elseif GLCT(r,c) == 16; SuitabilityFactor(r,c) = 0.25;
            elseif GLCT(r,c) == 17; SuitabilityFactor(r,c) = 0.25;
            elseif GLCT(r,c) == 18; SuitabilityFactor(r,c) = 0;
            elseif GLCT(r,c) == 19; SuitabilityFactor(r,c) = 0;
            end
            
            if abonAgLand(r,c) ==1
                SuitabilityFactor(r,c) = 0.9;
            end

        end
    end
end


%% Altmult 
altitudecutoff = 2000;
for r=1:nr
    for c=1:nc
        if Alt(r,c)<=altitudecutoff
            AltMult(r,c)=1;
        else
            AltMult(r,c)=0;
        end
    end
end

% figure(1);clf;imagesc(AltMult); colorbar

%% Biomult if zero than oke
for r=1:nr
    for c=1:nc
        if BIOres(r,c)==0
            BioMult(r,c)=1;
        else
            BioMult(r,c)=0;
        end
    end
end

% figure(1);clf;imagesc(BioMult); colorbar

%% Builtup
for r=1:nr
    for c=1:nc
        BuildupInv(r,c)=min(1, max(0, 1 - Buildup(r,c)));

        if BuildupInv(r,c) < 0.9
            BuildupInv(r,c) = 0;
        end
    end
end

% figure(1);clf;imagesc(BuildupInv); colorbar

%% Depth exclusion
for r=1:nr
    for c=1:nc
        if B(r,c)<-1000
            BatMult(r,c)=0;
        else
            BatMult(r,c)=1;
        end
    end
end

% figure(1);clf;imagesc(BatMult); colorbar

%% Distance to shore exclusion
for r=1:nr
    for c=1:nc
        if OffDis(r,c)>5
            DisMult(r,c)=0;
        else
            DisMult(r,c)=1;
        end
    end
end

% figure(1);clf;imagesc(DisMult); colorbar

%% ExclFactor onshore
for r=1:nr
    for c=1:nc
        ExclFactor(r,c) = SuitabilityFactor(r,c) * AltMult(r,c) * BioMult(r,c) * BuildupInv(r,c) ... %onshore
            * Ice(r,c) * Ship(r,c) * MP(r,c) * BatMult(r,c) * DisMult(r,c); %offshore
    end
end
% figure(2);clf;imagesc(ExclFactor); colormap(jet); colorbar

% Wind Exclusion factor
% txt=sprintf('ncols\t%d\nnrows\t%d\nxllcorner\t%d\nyllcorner\t%d\ncellsize\t%0.2f\nNODATA_value\t%d',720,360,-180,-90,0.5,-99);
% file = fullfile(root, sprintf('\\output\\Suitability_map.asc'));
% dlmwrite(file,txt,'');
% dlmwrite(file,ExclFactor,'-append','delimiter',' '); % 0-1 suitability factor

win=0;
a = EEZ==11;
cols = find(sum(a)>0);
rows = find(sum(a')>0);
fc = cols(1) - win;
lc = cols(end) + win;
fr = rows(1) - win;
lr = rows(end) + win;

b=ExclFactor;
b(find(~a))=-1;

bw = b(fr:lr,fc:lc);

% figure(2);clf;imagesc(bw);colorbar; axis image
