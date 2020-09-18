%% Economic potential
%Be aware that the onshore wind script is different from the offshore and solar scripts as it
%thinks from the perspective of the turbine, not from the perspective of the kW-installed

%%
PRef = 800;
beta = -0.3;
TurbineCostRef = 1050;
TurbShare = 0.8; % share of turbine cost of totoal installation cost
Interest = 0.1; %10
EconLifeTime = 20; %20 years
AnnualHours	= 8760; % Number of hours in a year
NewTransmissionCost_perKM = 1843000; % Ontario case study. Veileder is 2x cheaper
% NewTransmissionCost_perKM = 2390000; % Cost per km of building a new 230 kV double circuit transmission line (Source: Mason et al 2012)
AnnuityFactor = Interest/(1-(1+Interest)^-EconLifeTime);
AnnuityFactorTrans = Interest/(1-(1+Interest)^-60);
OnMShare = 0.03;
TransInvCostperkm = 160000;

%% Turbine cost
TurbineCost = ((Pnom / PRef) ^ beta) * TurbineCostRef;

%% Investment cost
InvestmentCost = (TurbineCost/TurbShare)*Pnom;

%% Annual investment cost
AnnualInvestmentCost = AnnuityFactor * InvestmentCost;

%% OM share
OnM = OnMShare * InvestmentCost;

%% Total annual cost
TotalAnnualCost = AnnualInvestmentCost + OnM;
TotalmonthlyCost = TotalAnnualCost/12;

%% Transmission cost
WindTransEnergy = 100000.0*8760;
WindTransLoss = 0.0005 * 0.05;
WindTransCostperkWhperkm = TransInvCostperkm * AnnuityFactorTrans /WindTransEnergy + WindTransLoss;

for r=1:nr
    for c=1:nc
        WindTransCost(r,c) = WindTransCostperkWhperkm * DistanceToLoad(r,c) + 0.01;
    end
end

% WindTransCost(WindTransCost>0.1)=0.1;
% figure(1);clf;imagesc(WindTransCost);axis image; colormap(jet);colorbar

%% Cost of Electricity $/kWh
for m=1:13
    for r=1:nr
        for c=1:nc
            if EloutTCorrect{m}(r,c) > EPS
                
                if setOnshorePL==1 %powerline setting distance cost
                    COE{m}(r,c) = WindTransCost(r,c) + (TotalAnnualCost*(1/12))/ EloutTCorrect{m}(r,c);
                else
                    COE{m}(r,c) = (TotalAnnualCost*(1/12))/ EloutTCorrect{m}(r,c);
                end
                
                if m==13
                    if setOnshorePL==1
                        COE{m}(r,c) = WindTransCost(r,c) + (TotalAnnualCost)/ EloutTCorrect{m}(r,c);
                    else
                        COE{m}(r,c) = (TotalAnnualCost)/ EloutTCorrect{m}(r,c);
                    end
                end
                
            else
                
                if setOnshorePL==1
                    COE{m}(r,c) = WindTransCost(r,c) + (TotalAnnualCost*(1/12))/ 1;
                else
                    COE{m}(r,c) = (TotalAnnualCost*(1/12))/ 1;
                end
                
                if m==13
                    if setOnshorePL==1
                        COE{m}(r,c) = WindTransCost(r,c) + (TotalAnnualCost)/ 1;
                    else
                        COE{m}(r,c) = (TotalAnnualCost)/ 1;
                    end
                end
                
            end
        end
    end
end

% COE{13}(COE{13}>0.2)=0.2;
% figure(1);clf;imagesc(COE{13});axis image; colormap(jet);colorbar

%% Cost of Electricity no cut $/kWh
for m=1:13
    for r=1:nr
        for c=1:nc
            if EloutTCorrNoCut{m}(r,c) > EPS
                
                if setOnshorePL==1 %powerline setting distance cost
                    COE_nocut{m}(r,c) = WindTransCost(r,c) + (TotalAnnualCost*(1/12))/ EloutTCorrNoCut{m}(r,c);
                else
                    COE_nocut{m}(r,c) = (TotalAnnualCost*(1/12))/ EloutTCorrNoCut{m}(r,c);
                end
                
                if m==13
                    if setOnshorePL==1
                        COE_nocut{m}(r,c) = (TotalAnnualCost)/ EloutTCorrNoCut{m}(r,c);
                    else
                        COE_nocut{m}(r,c) = (TotalAnnualCost)/ EloutTCorrNoCut{m}(r,c);
                    end
                end
                
            else
                
                if setOnshorePL==1
                    COE_nocut{m}(r,c) = WindTransCost(r,c) + (TotalAnnualCost*(1/12))/ 1;
                else
                    COE_nocut{m}(r,c) = (TotalAnnualCost*(1/12))/ 1;
                end
                
                if m==13
                    if setOnshorePL==1
                        COE_nocut{m}(r,c) = WindTransCost(r,c) + (TotalAnnualCost)/ 1;
                    else
                        COE_nocut{m}(r,c) = (TotalAnnualCost)/ 1;
                    end
                end
            end
        end
    end
end

% COE_nocut{13}(COE_nocut{13}>0.2)=0.2;
% figure(1);clf;imagesc(COE_nocut{13});axis image; colormap(jet);colorbar

%% Offshore calculation
% Regression formula from Offshore wind data, see Gernaat et al. 2014 Global long term cost developments of offshore wind
% Investment (€/kW) = B + a1 * DistanceFromShore + a2 * Depth + a3 * CumulativeCapacity
% Investment (€/kW) = 1326.23 + 7.72 * DistanceFromShore + 24.25 * Depth + 0.42 CumulativeCapacity

C1 = 1320; % Constant
B1 = 7.62; % Beta for Distance to shore (km)
B2 = 24.4; % Beta for Depth (m)
B3 = 0.42; % Beta for cumulative capacity (MW)
CumOffCap = 816; % MW Cumulative capacity (2006) consistent with the Junginger dataset
Eur2006toDollar2005 = 0.7838;

for r=1:nr
    for c=1:nc
        SpecCapCost(r,c) = C1 + B1*(50*OffDis(r,c)) + B2*(-1*B(r,c)) + B3*CumOffCap; % % EUR2006/kW
        SpecCapCostD(r,c) = SpecCapCost(r,c) * Eur2006toDollar2005; % $2005/kW
    end
end

% figure(4);clf;imagesc(SpecCapCostD);axis image; colormap(parula);colorbar; hold on

%% O&M
OnMShare = 0.15;

for r=1:nr
    for c=1:nc
        OnMo(r,c) = OnMShare * SpecCapCostD(r,c);
        TotalAnnualCost(r,c) = AnnuityFactor * SpecCapCostD(r,c) + OnMo(r,c); % Total Annual costs with O&M costs
    end
end

% figure(1);clf;imagesc(TotalAnnualCost);axis image; colormap(jet);colorbar

%% Electricity production of 1kW per year
for m=1:13
    for r=1:nr
        for c=1:nc
            ElecOutKW{m}(r,c) = Fulload{m}(r,c) * 1; % kWh produced if 1kW was installed
        end
    end
end

% figure(1);clf;imagesc(ElecOutKW{13});axis image; colormap(jet);colorbar
% figure(2);clf;imagesc(Fulload{13});axis image; colormap(jet);colorbar

%% Offshore COE
% Cost of Electricity (COE) is the total annual costs devided by the total electricity produced in a year
for m=1:13
    for r=1:nr
        for c=1:nc
            COEo{m}(r,c) = TotalAnnualCost(r,c)/ElecOutKW{m}(r,c); % $/kWh
        end
    end
end


% COEo{13}(COEo{13}>0.5)=0.5;
% figure(3);clf;
% ax1=subplot(1,2,1);imagesc(COE{13});axis image; colormap(jet);colorbar;title('COE Hoogwijk')
% ax2=subplot(1,2,2);imagesc(COEo{13});axis image; colormap(jet);colorbar;title('COE Gernaat')
% linkaxes([ax1 ax2])

%% Remove all the cells that have zero potential as defined in TechPot

Wind_LoadFac_prep = LoadFactor;
for m=1:13
    for r=1:nr
        for c=1:nc
            if TechPotCell{m}(r,c)==0;
                %onshore
                COE{m}(r,c)=NaN;
                COE_nocut{m}(r,c)=NaN;
                
                %Offshore
                COEo{m}(r,c)=NaN;
                
                if m==13
                    Wind_LoadFac_prep(r,c)=NaN;
                end
            end
        end
    end
end

TechPotCellCon = TechPotCell;
TechPotCellCoff = TechPotCell;
if setEconCeil==1
    for m=1:13
        %Onhore
        COE{m}(COE{m}>CostCeil)=NaN;
        COE_nocut{m}(COE_nocut{m}>CostCeil)=NaN;
        TechPotCellCon{m}(find(isnan(COE{m})))=0;
        
        %Offshore
        COEo{m}(COEo{m}>CostCeil)=NaN;
        TechPotCellCoff{m}(find(isnan(COEo{m})))=0;
        
        if m==13
            Wind_LoadFac_prep(COE{m}>CostCeil)=NaN;
        end
    end
    
    for m=1:13 %onshore
        for i=1:26
            RegTechPotWindEcon{m}(i) = sum(TechPotCellCon{m}(IRind{i})); % kWh / region / y
        end
        RegTechPotWindEcon{m}(27) =sum(RegTechPotWindEcon{m}(1:26));
    end
    
    for m=1:13 %onshore
        for i=1:numel(ISOGDP(:,1))
            CTechPotWindEcon{m}(i) = sum(TechPotCellCon{m}(Cind{i})); %kWh / country / y
        end
        dt = sum(CTechPotWindEcon{m}(1:end));
        CTechPotWindEcon{m}(end+1) =dt;
    end
    
    for m=1:13 %offshore
        for i=1:26
            RegTechPotWindEcono{m}(i) = sum(TechPotCellCoff{m}(IRindo{i})); % kWh / region / y
        end
        RegTechPotWindEcono{m}(27) =sum(RegTechPotWindEcono{m}(1:26));
    end
    
    
    for m=1:13 %offshore
        for i=1:numel(ISOGDP(:,1))
            CTechPotWindEcono{m}(i) = sum(TechPotCellCoff{m}(Cindo{i})); %kWh / country / y
        end
        dt = sum(CTechPotWindEcono{m}(1:end));
        CTechPotWindEcono{m}(end+1) =dt;
    end
end

%%
% COE{13}(COE{13}>0.5)=0.5;
% figure(1);clf;
% ax1=subplot(1,2,1);imagesc(COE{13});axis image; colorbar; title('COE ($/kWh)')
% ax2=subplot(1,2,2);imagesc(LoadFactor);axis image; colorbar; title('Load factor')
% linkaxes([ax1 ax2])

%% WindCost curve prep
clear SI
%Find cell per Image reg
for m=1:13
    for i=1:26
        
        %Onshore
        COE_IR{m}{i} = COE{m}(IRind{i});
        Pot_IR{m}{i} = TechPotCell{m}(IRind{i});
        
        %Offshore
        COEo_IR{m}{i} = COEo{m}(IRindo{i});
        Poto_IR{m}{i} = TechPotCell{m}(IRindo{i});
        
        if m==13
            LF_IR{i} = Wind_LoadFac_prep(IRind{i});  %onshore
            LFo_IR{i} = Wind_LoadFac_prep(IRindo{i});%offshore
        end
    end
end

%World
for m=1:13
    %onshore
    COE_IR{m}{27} = COE{m}(IRind{27});
    Pot_IR{m}{27} = TechPotCell{m}(IRind{27});
    
    %offshore
    COEo_IR{m}{27} = COEo{m}(IRindo{29});
    Poto_IR{m}{27} = TechPotCell{m}(IRindo{29});
    
    if m==13
        LF_IR{27} = Wind_LoadFac_prep(IRind{27});  %onshore
        LFo_IR{27} = Wind_LoadFac_prep(IRindo{29});%offshore
    end
end

% figure(1);clf;imagesc(COEo{13});

%% Sorting
for m=1:13
    for i=1:27
        %onshore
        [COE_IR_S{m}{i}, SI{m}{i}] = sort(COE_IR{m}{i});
        Pot_IR_S{m}{i} = Pot_IR{m}{i}(SI{m}{i});
        
        %offshore
        [COEo_IR_S{m}{i}, SIo{m}{i}] = sort(COEo_IR{m}{i});
        Poto_IR_S{m}{i} = Poto_IR{m}{i}(SIo{m}{i});
        
        if m==13
            LF_IR_S{i} = LF_IR{i}(SI{m}{i}); %onshore
            LFo_IR_S{i} = LFo_IR{i}(SIo{m}{i}); %offshore
        end
    end
end

% for m=1:13
%     [COE_IR_S{m}{27}, SI{m}{27}] = sort(COE_IR{m}{27});
%     [COEo_IR_S{m}{27}, SIo{m}{27}] = sort(COEo_IR{m}{27});
%     Pot_IR_S{m}{27} = Pot_IR{m}{27}(SI{m}{27});
%     Poto_IR_S{m}{27} = Pot_IR{m}{27}(SIo{m}{27});
%     if m==13
%         LF_IR_S{27} = LF_IR{27}(SI{m}{27});
%         LFo_IR_S{27} = LF_IR{27}(SIo{m}{27});
%     end
% end

%% Removing Nans Infs
for m=1:13
    for i=1:27
        %onshore NaNs
        COE_IR_S{m}{i} = COE_IR_S{m}{i}(~isnan(COE_IR_S{m}{i}));
        Pot_IR_S{m}{i} = Pot_IR_S{m}{i}(~isnan(COE_IR_S{m}{i}));
        LF_IR_S{i} = LF_IR_S{i}(~isnan(LF_IR_S{i}));
        
        %offshore
        COEo_IR_S{m}{i} = COEo_IR_S{m}{i}(~isnan(COEo_IR_S{m}{i}));
        Poto_IR_S{m}{i} = Poto_IR_S{m}{i}(~isnan(COEo_IR_S{m}{i}));
        LFo_IR_S{i} = LFo_IR_S{i}(~isnan(LFo_IR_S{i}));
        
        %onshore Infs
        COE_IR_S{m}{i} = COE_IR_S{m}{i}(find(~isinf(COE_IR_S{m}{i})));
        Pot_IR_S{m}{i} = Pot_IR_S{m}{i}(find(~isinf(COE_IR_S{m}{i})));
        LF_IR_S{i} = LF_IR_S{i}(find(~isinf(LF_IR_S{i})));
        
        %offshore
        COEo_IR_S{m}{i} = COEo_IR_S{m}{i}(find(~isinf(COEo_IR_S{m}{i})));
        Poto_IR_S{m}{i} = Poto_IR_S{m}{i}(find(~isinf(COEo_IR_S{m}{i})));
        LFo_IR_S{i} = LFo_IR_S{i}(find(~isinf(LFo_IR_S{i})));
        
    end
end
