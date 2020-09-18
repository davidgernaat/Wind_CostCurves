%%
fulloadfactor1 = 565;
fulloadfactor2 = 1745;
MinFulload = 515;
FullloadMax = 4200;
MaxkWPerkm2 = 4000;
OpEff = 0.92; %in 2010, goes to 0.99 in 2100
ArrayEff = 0.914; %in 2010, goes to 0.97 in 2100
EPS = 0.000001;

%% Loadfactor with cut
% in kWh/kW fulload is equal to the power * power curve* weibull depending on K factor. using a linear relation (see documentation)
for m=1:13
    Fulload{m} = zeros(size(V{13}));
    for r=1:nr
        for c=1:nc
            
            if (fulloadfactor1 * Vcut_Hub{m}(r,c) - fulloadfactor2) < MinFulload
                Fulload{m}(r,c) = NaN;
            else
                Fulload{m}(r,c) = fulloadfactor1 * Vcut_Hub{m}(r,c) - fulloadfactor2;
            end
            
            Fulload{m}(r,c) = min(FullloadMax,Fulload{m}(r,c));
            
            if m<13 %Not for annual
                Fulload{m}(r,c) = (1/12) * Fulload{m}(r,c);
            end
        end
    end
end

% figure(1);clf;imagesc(Fulload{13}); colormap(jet); colorbar; axis image

%% Loadfactor without cut
for m=1:13
    for r=1:nr
        for c=1:nc
            
            FulloadNoCut{m}(r,c) = (1/12) * min(8760,max(0,fulloadfactor1 * VHub{m}(r,c) - fulloadfactor2));
            
            if m==13 %Not for annual
                FulloadNoCut{m}(r,c) = min(8760,max(0,fulloadfactor1 * VHub{m}(r,c) - fulloadfactor2));
            end
            
        end
    end
end

% figure(1);clf;imagesc(FulloadNoCut{13}); colormap(jet); colorbar; axis image

%% Define loadfactor
for r=1:nr
    for c=1:nc
        LoadFactor(r,c) = Fulload{13}(r,c)/8760;
        LoadFactorNoCut(r,c) = FulloadNoCut{13}(r,c)/8760;
    end
end

% figure(1);clf;
% subplot(2,1,1); imagesc(LoadFactor{13}); colormap(jet); colorbar; axis image
% subplot(2,1,2); imagesc(LoadFactorNoCut{13}); colormap(jet); colorbar; axis image

%% Turbines per cell
% To be clear, this is for all area, so also, the once with vavg < 4.0 m/S, therefore is corrected with the technical potential

TurbSpace = 1/MaxkWPerkm2 * Pnom;

for r=1:nr
    for c=1:nc
        TurbPerCell(r,c)	= (Area(r,c) * ExclFactor(r,c)) / TurbSpace;
        TurbPerCellTheo(r,c)	= Area(r,c) / TurbSpace;
    end
end

% figure(1);clf;imagesc(Area); colormap(jet); colorbar; axis image

%% Eelectricity output of a turbine using the simpel approach, kWh
for m=1:13
    for r=1:nr
        for c=1:nc
            EloutT{m}(r,c) = Fulload{m}(r,c) * Pnom;
        end
    end
end

for m=1:13 %onshore
    for i=1:26
        ElouttRegion{m}(i) = sum(EloutT{m}(IRind{i})); % kWh / region / y
    end
    ElouttRegion{m}(27) =sum(ElouttRegion{m}(1:26));
end

for m=1:13 %offshore
    for i=1:26
        ElouttRegiono{m}(i) = sum(EloutT{m}(IRindo{i})); % kWh / region / y
    end
    ElouttRegiono{m}(27) =sum(ElouttRegiono{m}(1:26));
end

%% Eelectricity output of a turbine using the simpel approach no cut, kWh
for m=1:13
    for r=1:nr
        for c=1:nc
            EloutTNoCut{m}(r,c) = FulloadNoCut{m}(r,c) * Pnom;
        end
    end
end

%% kWh per turbine
for m=1:13
    for r=1:nr
        for c=1:nc
            EloutTCorrect{m}(r,c) = EloutT{m}(r,c) * OpEff * ArrayEff;
        end
    end
end

%% kWh per turbine no cut
for m=1:13
    for r=1:nr
        for c=1:nc
            EloutTCorrNoCut{m}(r,c) = EloutTNoCut{m}(r,c) * OpEff * ArrayEff;
        end
    end
end

%% Technical potential kWh
for m=1:12
    for r=1:nr
        for c=1:nc
            TechPotCell{m}(r,c) = ExclFactor(r,c) * EloutTCorrect{m}(r,c) * TurbPerCell(r,c) * AreaVcut{m}(r,c);  % kWh / cell / y
        end
    end
end

%Total
for r=1:nr
    for c=1:nc
        for m=1:12
            dt(m) = TechPotCell{m}(r,c);
        end
        TechPotCell{13}(r,c) = sum(dt); % kWh / cell / y
    end
end

for m=1:13 %onshore
    for i=1:26
        RegTechPotWind{m}(i) = sum(TechPotCell{m}(IRind{i})); % kWh / region / y
    end
    RegTechPotWind{m}(27) =sum(RegTechPotWind{m}(1:26));
end

for m=1:13 %onshore
    for i=1:numel(ISOGDP(:,1))
        Cind{i} = find(GISO(:)==ISOGDP(i,1));
        CTechPotWind{m}(i) = sum(TechPotCell{m}(Cind{i})); %kWh / country / y
    end
    dt = sum(CTechPotWind{m}(1:end));
    CTechPotWind{m}(end+1) =dt;
end

for m=1:13 %offshore
    for i=1:26
        RegTechPotWindo{m}(i) = sum(TechPotCell{m}(IRindo{i})); % kWh / region / y
    end
    RegTechPotWindo{m}(27) =sum(RegTechPotWindo{m}(1:26));
end

for m=1:13 %offshore
    for i=1:numel(ISOGDP(:,1))
        Cindo{i} = find(CEEZ(:)==ISOGDP(i,1));
        CTechPotWindo{m}(i) = sum(TechPotCell{m}(Cindo{i})); %kWh / country / y
    end
    dt = sum(CTechPotWindo{m}(1:end));
    CTechPotWindo{m}(end+1) =dt;
end

% figure(1);clf;imagesc(TechPotCell{13});axis image; colorbar; colormap(cmap)

%% Technical potential no cut kWh
for m=1:12
    for r=1:nr
        for c=1:nc
            TechPotCellNoCut{m}(r,c) = ExclFactor(r,c) * EloutTCorrNoCut{m}(r,c) * TurbPerCell(r,c);
        end
    end
end

%Total
for r=1:nr
    for c=1:nc
        for m=1:12
            dt(m) = TechPotCellNoCut{m}(r,c);
        end
        TechPotCellNoCut{13}(r,c) = sum(dt);
    end
end

for m=1:13 %onshore
    for i=1:26
        RegTechPotWindNoCut{m}(i) = sum(TechPotCellNoCut{m}(IRind{i})); % kWh / region / y
    end
    RegTechPotWindNoCut{m}(27) =sum(RegTechPotWindNoCut{m}(1:26));
end

for m=1:13 %onshore
    for i=1:numel(ISOGDP(:,1))
        CTechPotWindNoCut{m}(i) = sum(TechPotCellNoCut{m}(Cind{i})); %kWh / country / y
    end
    dt = sum(CTechPotWindNoCut{m}(1:end));
    CTechPotWindNoCut{m}(end+1) =dt;
end

for m=1:13 %offshore
    for i=1:26
        RegTechPotWindNoCuto{m}(i) = sum(TechPotCellNoCut{m}(IRindo{i})); % kWh / region / y
    end
    RegTechPotWindNoCuto{m}(27) =sum(RegTechPotWindNoCuto{m}(1:26));
end

for m=1:13 %offshore
    for i=1:numel(ISOGDP(:,1))
        CTechPotWindNoCuto{m}(i) = sum(TechPotCellNoCut{m}(Cindo{i})); %kWh / country / y
    end
    dt = sum(CTechPotWindNoCuto{m}(1:end));
    CTechPotWindNoCuto{m}(end+1) =dt;
end

% figure(1);clf;imagesc(COE{13});axis image; colorbar

%% Technical potential theoretical KWh
for m=1:13
    for r=1:nr
        for c=1:nc
            TechPotCellTheo{m}(r,c) = EloutTCorrNoCut{m}(r,c) * TurbPerCellTheo(r,c); % kWh / cell / y
        end
    end
end

for m=1:13 %onshore
    for i=1:26
        RegTechPotWindTheo{m}(i) = sum(TechPotCellTheo{m}(IRind{i})); % kWh / region / y
    end
    RegTechPotWindTheo{m}(27) =sum(RegTechPotWindTheo{m}(1:26));
end

for m=1:13 %onshore
    for i=1:numel(ISOGDP(:,1))
        CTechPotWindTheo{m}(i) = sum(TechPotCellTheo{m}(Cind{i})); %kWh / country / y
    end
    dt = sum(CTechPotWindTheo{m}(1:end));
    CTechPotWindTheo{m}(end+1) =dt;
end

for m=1:13 %offshore
    for i=1:26
        RegTechPotWindTheoo{m}(i) = sum(TechPotCellTheo{m}(IRindo{i})); % kWh / region / y
    end
    RegTechPotWindTheoo{m}(27) =sum(RegTechPotWindTheoo{m}(1:26));
end

for m=1:13 %offshore
    for i=1:numel(ISOGDP(:,1))
        CTechPotWindTheoo{m}(i) = sum(TechPotCellTheo{m}(Cindo{i})); %kWh / country / y
    end
    dt = sum(CTechPotWindTheoo{m}(1:end));
    CTechPotWindTheoo{m}(end+1) =dt;
end

%% Technical potential Geographical KWh
for m=1:13
    for r=1:nr
        for c=1:nc
            TechPotCellGeo{m}(r,c) = ExclFactor(r,c) * EloutTCorrNoCut{m}(r,c) * TurbPerCellTheo(r,c); % kWh / cell / y
        end
    end
end

for m=1:13 %onshore
    for i=1:26
        RegTechPotWindGeo{m}(i) = sum(TechPotCellGeo{m}(IRind{i})); % kWh / region / y
    end
    RegTechPotWindGeo{m}(27) =sum(RegTechPotWindGeo{m}(1:26));
end

for m=1:13 %onshore
    for i=1:numel(ISOGDP(:,1))
        CTechPotWindGeo{m}(i) = sum(TechPotCellGeo{m}(Cind{i})); %kWh / country / y
    end
    dt = sum(CTechPotWindGeo{m}(1:end));
    CTechPotWindGeo{m}(end+1) =dt;
end

for m=1:13 %offshore
    for i=1:26
        RegTechPotWindGeoo{m}(i) = sum(TechPotCellGeo{m}(IRindo{i})); % kWh / region / y
    end
    RegTechPotWindGeoo{m}(27) =sum(RegTechPotWindGeoo{m}(1:26));
end

for m=1:13 %offshore
    for i=1:numel(ISOGDP(:,1))
        CTechPotWindGeoo{m}(i) = sum(TechPotCellGeo{m}(Cindo{i})); %kWh / country / y
    end
    dt = sum(CTechPotWindGeoo{m}(1:end));
    CTechPotWindGeoo{m}(end+1) =dt;
end

%% Total active area
for r=1:nr
    for c=1:nc
        if TechPotCell{13}(r,c) > EPS
            AreaAct(r,c) = Area(r,c) * SuitabilityFactor(r,c);
        else
            AreaAct(r,c) = 0;
        end
    end
end

for i=1:26
    TotArea(i) = sum(Area(IRind{i})); % m2 per region active onshore
end
TotArea(27) =sum(TotArea(1:26));

for i=1:26
    TotAreao(i) = sum(Area(IRindo{i})); % m2 per region active offshore
end
TotAreao(27) =sum(TotAreao(1:26));

for i=1:26
    TotAreaAct(i) = sum(AreaAct(IRind{i})); % m2 per region onshore
end
TotAreaAct(27) =sum(TotAreaAct(1:26));

for i=1:26
    TotAreaActo(i) = sum(AreaAct(IRindo{i})); % m2 per region offshore
end
TotAreaActo(27) =sum(TotAreaActo(1:26));

%%
% fname = sprintf('%s\\Wind_CC\\TechpotCellMonthTotal.out', root);
% [TechpotCellMonthTotal_data] = read_mym(fname);
% 
% TechpotCellMonthTotal = zeros(360,720);
% [nr nc] = size(TechpotCellMonthTotal);
% 
% i=0;
% for r=1:nr
%     for c=1:nc
%         if imagemask(r,c)==0; continue; end;
%         i=i+1;
%         TechpotCellMonthTotal(r,c)=TechpotCellMonthTotal_data(i);
%     end
% end
% 
% figure(1);clf;
% ax1=subplot(1,2,1);imagesc(TechpotCellMonthTotal);axis image; title('Hoogwijk')
% ax2=subplot(1,2,2);imagesc(TechPotCell{13});axis image; title('Mine')
% linkaxes([ax1 ax2])

%%
% clear dtm dth
% i=0;
% for r=1:nr
%     for c=1:nc
%         if TechPotCell{13}(r,c) > EPS
%             dtm(r,c)=1;
%         else
%             dtm(r,c)=0;
%         end
%         if TechpotCellMonthTotal(r,c) > EPS
%             dth(r,c)=1;
%         else
%             dth(r,c)=0;
%         end
%     end
% end
% for r=1:nr
%     for c=1:nc
%         j=j+1;
%         if imagemask(r,c)==0; 
%             TechpotCellMonthTotal(r,c)=-1000000000;
%             TechPotCell{13}(r,c)=-10000000000;
%         end
%     end
% end
% cmap=jet(256);
% cmap(1,:)=[1 1 1];
% cmap(2,:)=[0.9 0.9 0.9];
% figure(1);clf;
% ax1=subplot(2,2,1);imagesc(TechpotCellMonthTotal);axis image; title('Hoogwijk');colormap(cmap)
% ax2=subplot(2,2,2);imagesc(TechPotCell{13});axis image; title('Mine');colormap(cmap)
% ax3=subplot(2,2,3);imagesc(dth);axis image; title('Hoogwijk')
% ax4=subplot(2,2,4);imagesc(dtm);axis image; title('Mine')
% linkaxes([ax1 ax2 ax3 ax4])
% sum(TechpotCellMonthTotal(:))*1e-12
% sum(TechPotCell{13}(:))*1e-12

%% Which IMAGE cell is which matlab cell
% [nr nc] = size(imagemask);
% mc = zeros([nr nc]);
% i=0; j=0; r=1;
% for r=1:nr
%     for c=1:nc
%         j=j+1;
%         if imagemask(r,c)==1; 
%             i=i+1; 
%             ic(j)=i;
%             mc(r,c)=i;
%         else
%             ic(j)=0;
%         end
%         
%     end
% end
% a=find(ic);
% 
% % IMAGE2Matlab
% % Whats image cell number? 
% x=70;
% rr = ceil(a(x)/nc)
% cc = a(x) - (nc*(rr-1))
% 
% % figure(2);clf;imagesc(imagemask);colormap(flipud(gray)); hold on
% % plot(cc,rr,'.r','markersize',20)
% % hold off
% 
% % Matlab2IMAGE
% % Whatare the row and columm numbers?
% rr = 112;
% cc = 180;
% m=1;
% x = mc(rr,cc)

%% Offshore wind potential Germany
% figure(1);clf;imagesc(EEZ);colorbar; axis image
% a = CEEZ==224;
% b=TechPotCell{13};
% b(find(~a))=0;
% figure(1);clf;imagesc(b);colorbar; axis image
% 
% fprintf('Offshore wind potential: %0.2f GWh\n',sum(bw(:))*1e-6)
% fprintf('Current generation Germany: %0.2f GWh\n',19341)

%% Onshore wind potential Germany
% figure(1);clf;imagesc(GISO);colorbar; axis image
% a = GISO==276;
% b=TechPotCell{13};
% b(find(~a))=0;
% figure(1);clf;imagesc(b);colorbar; axis image
% 
% fprintf('Onshore wind potential: %0.2f TWh\n',sum(bw(:))*1e-9)
% fprintf('Current generation Germany: %0.2f TWh\n',92.249)

%%
% figure(1);clf;imagesc(EEZ);colorbar; axis image
% win=50;
% a = CEEZ==224;
% cols = find(sum(a)>0);
% rows = find(sum(a')>0);
% fc = cols(1) - win;
% lc = cols(end) + win;
% fr = rows(1) - win;
% lr = rows(end) + win;
% 
% b=TechPotCell{13};
% b(find(~a))=0;
% 
% bw = b(fr:lr,fc:lc);
% 
% fprintf('Total potential: %0.2f GWh\n',sum(bw(:))*1e-6)
% fprintf('Current generation Germany: %0.2f GWh\n',19341)
% 
% figure(1);clf;imagesc(bw);colorbar; axis image
