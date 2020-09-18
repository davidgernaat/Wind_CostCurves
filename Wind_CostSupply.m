% calcCostCurve  (region,NR27,COEMonthAverage,TechPotCellMonthTotal,NPC,PotentialMonthAverage,Pointer);

% calcCostCurve  (COUNTRY,NR224,COEMonthAverage,TechPotCellMonthTotal,NPC,PotMonthAverageCountry,PointerCountry);

%% Transform into TIMER curves

for j=1:27
    
    %Wind onshore
    Onshore_Pot_cumsum{j} = cumsum(Pot_IR_S{13}{j}); %kWh
    
    for i=1:numel(Onshore_Pot_cumsum{j})
        Onshore_Pot_Indexed{j}(i) = Onshore_Pot_cumsum{j}(i)/Onshore_Pot_cumsum{j}(end); 
    end
    
    %Wind offshore
    Offshore_Pot_cumsum{j} = cumsum(Poto_IR_S{13}{j}); %kWh
    
    for i=1:numel(Offshore_Pot_cumsum{j})
        Offshore_Pot_Indexed{j}(i) = Offshore_Pot_cumsum{j}(i)/Offshore_Pot_cumsum{j}(end); 
    end
    
end

%% Find min idx
val = linspace(0,1,100);

for j=1:27
    for i=1:numel(val)
        
        %Wind onshore
        clear tmp
        if isempty(Onshore_Pot_Indexed{j})==1; continue; end
        tmp = abs((Onshore_Pot_Indexed{j}./val(i))-1);
        [~, idxminOnshore{j}(i)] = min(tmp); %index of closest value
        
        %Wind offshore
        clear tmp
        if isempty(Offshore_Pot_Indexed{j})==1; continue; end
        tmp = abs((Offshore_Pot_Indexed{j}./val(i))-1);
        [~, idxminOffshore{j}(i)] = min(tmp); %index of closest value
        
    end
end

%% Put in TIMER matrix

a = [0	0.01	0.02	0.03	0.04	0.05	0.06	0.07	0.08	0.09	0.1	0.11	0.12	0.13	0.14	0.15	0.16	0.17	0.18	0.19	0.2	0.21	0.22	0.23	0.24	0.25	0.26	0.27	0.28	0.29	0.3	0.31	0.32	0.33	0.34	0.35	0.36	0.37	0.38	0.39	0.4	0.41	0.42	0.43	0.44	0.45	0.46	0.47	0.48	0.49	0.5	0.51	0.52	0.53	0.54	0.55	0.56	0.57	0.58	0.59	0.6	0.61	0.62	0.63	0.64	0.65	0.66	0.67	0.68	0.69	0.7	0.71	0.72	0.73	0.74	0.75	0.76	0.77	0.78	0.79	0.8	0.81	0.82	0.83	0.84	0.85	0.86	0.87	0.88	0.89	0.9	0.91	0.92	0.93	0.94	0.95	0.96	0.97	0.98	0.99	1];

CostCurveSmthOnshore = zeros(100,28);
CostCurveSmthOffshore = zeros(100,28);

LFCurveSmthOnshore = zeros(100,28);
LFCurveSmthOffshore = zeros(100,28);

for i=1:numel(a)-1
    CostCurveSmthOnshore(i,1) = a(i);
    CostCurveSmthOffshore(i,1) = a(i);

    LFCurveSmthOnshore(i,1) = a(i);
    LFCurveSmthOffshore(i,1) = a(i);
end

for i=1:27
    COEISort_Indexed_Onshore{i} = COE_IR_S{13}{i}(idxminOnshore{i});
    COEISort_Indexed_Offshore{i} = COEo_IR_S{13}{i}(idxminOffshore{i});
    
    LFISort_Indexed_Onshore{i} = LF_IR_S{i}(idxminOnshore{i});
    LFISort_Indexed_Offshore{i} = LFo_IR_S{i}(idxminOffshore{i});
end

for i=1:27
    % Cost curves
    % Onshore
    if isempty(COEISort_Indexed_Onshore{i})==1; continue; end
    CostCurveSmthOnshore(:,i+1) = COEISort_Indexed_Onshore{i};
    
    % Offshore
    if isempty(COEISort_Indexed_Offshore{i})==1; continue; end
    CostCurveSmthOffshore(:,i+1) = COEISort_Indexed_Offshore{i};
    
    %Load curves 
    % Onshore
    if isempty(LFISort_Indexed_Onshore{i})==1; continue; end
    LFCurveSmthOnshore(:,i+1) = LFISort_Indexed_Onshore{i};
    
    % Offshore
    if isempty(LFISort_Indexed_Offshore{i})==1; continue; end
    LFCurveSmthOffshore(:,i+1) = LFISort_Indexed_Offshore{i};
end

for i=1:27
    % Cost curves
    % Onshore
    if numel(COEISort_Indexed_Onshore{i})<1
        CostCurveSmthOnshore(101,1) = 1;
        CostCurveSmthOnshore(101,i+1) = 1;
    else
        CostCurveSmthOnshore(101,1) = 1;
        CostCurveSmthOnshore(101,i+1) = COEISort_Indexed_Onshore{i}(end);
    end

    % Offshore
    if numel(COEISort_Indexed_Offshore{i})<1
        CostCurveSmthOffshore(101,1) = 1;
        CostCurveSmthOffshore(101,i+1) = 1;
    else
        CostCurveSmthOffshore(101,1) = 1;
        CostCurveSmthOffshore(101,i+1) = COEISort_Indexed_Offshore{i}(end);
    end
    
    % Load curves
    % Onshore
    if numel(COEISort_Indexed_Onshore{i})<1
        LFCurveSmthOnshore(101,1) = 1;
        LFCurveSmthOnshore(101,i+1) = 0;
    else
        LFCurveSmthOnshore(101,1) = 1;
        LFCurveSmthOnshore(101,i+1) = LFISort_Indexed_Onshore{i}(end);
    end
    
    % Offshore
    if numel(COEISort_Indexed_Offshore{i})<1
        LFCurveSmthOffshore(101,1) = 1;
        LFCurveSmthOffshore(101,i+1) = 0;
    else
        LFCurveSmthOffshore(101,1) = 1;
        LFCurveSmthOffshore(101,i+1) = LFISort_Indexed_Offshore{i}(end);
    end
end

%MaxProd
MaxProdOnshore = zeros(1,27);
MaxProdOffshore = zeros(1,27);

for i=1:27
    if numel(Onshore_Pot_cumsum{i})==0
        MaxProdOnshore(:,i) = 0;
    else
        MaxProdOnshore(:,i) = Onshore_Pot_cumsum{i}(end) * 0.0036; %GJ
    end
    
    if numel(Offshore_Pot_cumsum{i})==0
        MaxProdOffshore(:,i) = 0;
    else
        MaxProdOffshore(:,i) = Offshore_Pot_cumsum{i}(end) * 0.0036; %GJ
    end
end

%% Show
% Costcurves

% h=figure(1);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(CostCurveSmthOnshore(1:end-3,i+1));
%     MaxPotential(i) = max(MaxProdOnshore(i))/0.0036;
%     maxpot{i} = sprintf('On %s (%0.2f PWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('$2010/kWh','FontSize',8)
%     xlabel('Cumulative energy supply (PWh)','FontSize',8)
% end
% 
% %Loadfactor
% 
% h=figure(2);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(LFCurveSmthOnshore(:,i+1));
%     MaxPotential(i) = max(MaxProdOnshore(i))/0.0036;
%     maxpot{i} = sprintf('On %s (%0.2f PWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('Load factor','FontSize',8)
%     xlabel('Cumulative energy supply (PWh)','FontSize',8)
% end
% 
% h=figure(3);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(CostCurveSmthOffshore(1:end-3,i+1));
%     MaxPotential(i) = max(MaxProdOffshore(i))/0.0036;
%     maxpot{i} = sprintf('Off %s (%0.2f PWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('$2010/kWh','FontSize',8)
%     xlabel('Cumulative energy supply (PWh)','FontSize',8)
% end
% 
% %Loadfactor
% 
% h=figure(4);clf;
% for i=1:27
%     subplot(5,6,i)
%     stairs(LFCurveSmthOffshore(:,i+1));
%     MaxPotential(i) = max(MaxProdOffshore(i))/0.0036;
%     maxpot{i} = sprintf('Off %s (%0.2f PWh)\n',IMAGER{i},MaxPotential(i)*1e-12);
%     title(maxpot{i}, 'FontWeight','bold','FontSize',10)
%     ylabel('Load factor','FontSize',8)
%     xlabel('Cumulative energy supply (PWh)','FontSize',8)
% end

%% Save
if setoutput==true;

    pathname = fileparts(root);
    
    % CostCurves
    
    % Onshore
    file = fullfile(pathname, sprintf('\\Wind_CC\\output\\CostCurveSmthOnshore.dat'));
    txt=sprintf('real main.OnshoreCostCurveSmth[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,CostCurveSmthOnshore,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % MaxProd
    
    %Onshore
    file = fullfile(pathname, sprintf('\\Wind_CC\\output\\MaxProdOnshore.dat'));
    txt=sprintf('real main.Onshore_TechPotRegion[27] = ');
    dlmwrite(file,txt,'');
    dlmwrite(file,MaxProdOnshore,'-append','delimiter','\t');
    dlmwrite(file,';','-append');
    
    % Load factor decline curves
    
    %Onshore
    file = fullfile(pathname, sprintf('\\Wind_CC\\output\\LoadCurveSmthOnshore.dat'));
    txt=sprintf('real main.OnshoreLoadCurveSmth[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,LFCurveSmthOnshore,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
        % CostCurves
    
    % Offshore
    file = fullfile(pathname, sprintf('\\Wind_CC\\output\\CostCurveSmthOffshore.dat'));
    txt=sprintf('real main.OffshoreCostCurveSmth[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,CostCurveSmthOffshore,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
    % MaxProd
    
    %Offshore
    file = fullfile(pathname, sprintf('\\Wind_CC\\output\\MaxProdOffshore.dat'));
    txt=sprintf('real main.Offshore_TechPotRegion[27] = ');
    dlmwrite(file,txt,'');
    dlmwrite(file,MaxProdOffshore,'-append','delimiter','\t');
    dlmwrite(file,';','-append');
    
    % Load factor decline curves
    
    %Offshore
    file = fullfile(pathname, sprintf('\\Wind_CC\\output\\LoadCurveSmthOffshore.dat'));
    txt=sprintf('real main.OffshoreLoadCurveSmth[27](oo) = [');
    dlmwrite(file,txt,'');
    dlmwrite(file,LFCurveSmthOffshore,'-append','delimiter','\t');
    dlmwrite(file,'];','-append');
    
end