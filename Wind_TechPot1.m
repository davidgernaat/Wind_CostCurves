%% Technical potential CSP PV and PVres

Pnom = 1420; %1420 in 2000, in 2100: 3000
PHHfactor = 0.28;		% to assess hubheight based on nominal power
h10 = 10;               % height at 10 m

%% Hubheight
Hubheight = 10 * Pnom^PHHfactor;

%% Which wind file
if windfile==2
    for m=1:13
        wspr{m} = resizem(wsp{m},2);
        pct10mr{m} = resizem(pct10m{m},2);
    end
    
    if hubheightcorr==1
        [nr nc]=size(wspr{13});
        for m=1:13
            for r=1:nr
                for c=1:nc
                    wspr{m}(r,c) = wspr{m}(r,c)* (1-(-1*pct10mr{m}(r,c)*1e-2));
                end
            end
        end
    end
    
    V=wspr;
end

% figure(2);clf;
% ax1=subplot(1,2,1);imagesc(V{13});axis image; title('Hoogwijk')
% ax2=subplot(1,2,2);imagesc(wspr2{13});axis image; title('Nasa')
% linkaxes([ax1 ax2])
% figure(1);clf;imagesc(V{13});axis image;colorbar

%% Which GLCT file
if glctfile==1
    GLCT = glctwindat;
end

%% Roughnesslenght corrector
if glctfile==2
    RoughL = zeros(size(GLCT));
    for r=1:nr
        for c=1:nc
            if GLCT(r,c) == 1; RoughL(r,c) = 0.25;
            elseif GLCT(r,c) == 2; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 3; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 4; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 5; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 6; RoughL(r,c) = 0.25;
            elseif GLCT(r,c) == 7; RoughL(r,c) = 0.25;
            elseif GLCT(r,c) == 8; RoughL(r,c) = 0.25;
            elseif GLCT(r,c) == 9; RoughL(r,c) = 0.25;
            elseif GLCT(r,c) == 10; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 11; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 12; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 13; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 14; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 15; RoughL(r,c) = 0.03;
            elseif GLCT(r,c) == 16; RoughL(r,c) = 0.005;
            elseif GLCT(r,c) == 17; RoughL(r,c) = 0.1;
            elseif GLCT(r,c) == 18; RoughL(r,c) = 0.25;
            elseif GLCT(r,c) == 19; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 20; RoughL(r,c) = 1;
            end
        end
    end
else
    RoughL = zeros(size(GLCT));
    for r=1:nr
        for c=1:nc
            if GLCT(r,c) == 1; RoughL(r,c) = 0.25;
            elseif GLCT(r,c) == 2; RoughL(r,c) = 0.03;
            elseif GLCT(r,c) == 3; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 4; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 5; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 6; RoughL(r,c) = 0.005;
            elseif GLCT(r,c) == 7; RoughL(r,c) = 0.25;
            elseif GLCT(r,c) == 8; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 9; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 10; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 11; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 12; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 13; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 14; RoughL(r,c) = 0.1;
            elseif GLCT(r,c) == 15; RoughL(r,c) = 0.005;
            elseif GLCT(r,c) == 16; RoughL(r,c) = 0.1;
            elseif GLCT(r,c) == 17; RoughL(r,c) = 0.3;
            elseif GLCT(r,c) == 18; RoughL(r,c) = 1;
            elseif GLCT(r,c) == 19; RoughL(r,c) = 1;
            end
        end
    end
end

%% Set sea roughness factor to very low
for r=1:nr
    for c=1:nc
        if RoughL(r,c) == 0; RoughL(r,c) = 0.0002;
        end
    end
end
% figure(1);clf;imagesc(RoughL);

%%
for r=1:nr
    for c=1:nc
        Fac1(r,c) = log(Hubheight / max(0.00005, RoughL(r,c)));
        Fac2(r,c) = log(h10 / max(0.00005, RoughL(r,c)));
        Fac1100(r,c) = log(100 / max(0.005, RoughL(r,c)));
    end
end

%% V cutoff
for m=1:13
    for r=1:nr
        for c=1:nc
            if V{m}(r,c) <= 4
                Vcut{m}(r,c) = 0;
            else
                Vcut{m}(r,c) = V{m}(r,c);
            end
        end
    end
end

%% V cutoff at hub 
for m=1:13
    for r=1:nr
        for c=1:nc
            Vcut_Hub{m}(r,c) = Vcut{m}(r,c) * (Fac1(r,c)/Fac2(r,c));
        end
    end
end

%% V at hub
for m=1:13
    for r=1:nr
        for c=1:nc
            VHub{m}(r,c) = V{m}(r,c) * (Fac1(r,c)/Fac2(r,c));
        end
    end
end

%% Area V cutoff
for m=1:13
    for r=1:nr
        for c=1:nc
            if V{m}(r,c) <= 4
                AreaVcut{m}(r,c) = 0;
            else
                AreaVcut{m}(r,c) = 1;
            end
        end
    end
end

%% Which cells which regions onshore
for i=1:26
    IRind{i} = find(IRegion(:)==i);
end
IRind{27} = find(IRegion>0); %All onshore cells

%% Which cells which regions offshore
for i=1:28
    IRindo{i} = find(EEZ(:)==i);
end
IRindo{29} = find(EEZ>0 & EEZ<28); %All onshore cells
