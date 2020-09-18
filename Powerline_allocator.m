function DisCost = Powerline_allocator(D,PW,multifactor)

% D = 180;
% PW = 3658; 

%% Powerline allocator
% Determines which powerlines to build. Two sources for cost: (1) veileder
% and (2) Ontario case study. The second has 2x higher cost.

%% Log plot based on Veileder p143
% D=1:1000;
% P10 = max(1,min(4,12*D.^-1));
% P20 = max(1,min(8,48*D.^-1));
% P60 = max(1,min(22,520*D.^-1));
% P132 = max(1,min(125,2775.5*D.^-1));
% P300 = max(1,min(360,18000*D.^-1));
%
% figure(1);clf;
% h2 = plot(log10(D),log10(P20),'r');
% hold on
% h3 = plot(log10(D),log10(P60),'k');
% h4 = plot(log10(D),log10(P132),'g');
% h5 = plot(log10(D),log10(P300),'m');
% plot(log10(1),log10(10),'.r','markersize',20)
% hold off

%% Normal plot based on Veileder p143
% figure(2);clf;
% h1 = area(D,P10);
% % axis([0,1000,0,400])
% hold on
% h2 = plot(D,P20,'r');
% h3 = plot(D,P60,'k');
% h4 = plot(D,P132,'g');
% h5 = plot(D,P300,'m');
% plot(100,150,'.r','markersize',20)
% hold off

%%
%Check
% D=[1 10 50 100 500]; km
% PW=[10 10 10 10 10]; MW

PWx=1;
if PW >500;
    PWx = floor(PW/500);
end
if PWx==0; PWx=1; end

if isnan(PW);
    DisCost=NaN;
    return; 
end

PowerlineSelect=zeros(4,1);
P(1) = max(1,min(8,48*D.^-1)); %20kV
P(2) = max(1,min(22,520*D.^-1)); %60kV
P(3) = max(1,min(125,2775.5*D.^-1)); %132kV
P(4) = max(1,min(360,18000*D.^-1)); %300kV

for i=1:4
    if P(i)>PW ;
        PowerlineSelect(i)=1;
    end
end

PowerlineType = find(PowerlineSelect);
if isempty(PowerlineType); PowerlineType=5; end; %Double circuit lines

%% Cost Veileder
% PlineC(1) = 0.75*1e6; % NOK/km 24kV Simplex steel powerline 100 FeAl Veileder p139
% PlineC(2) = 1.3*1e6;  % NOK/km 72.5kV Simplex steel powerline 180 FeAl Veileder p140
% PlineC(3) = 1.7*1e6;  % NOK/km 145kV Simplex steel powerline 180 FeAl Veileder p141
% PlineC(4) = 3.15*1e6; % NOK/km 300kV Simplex steel powerline cracle Veileder p142
% PlineC(5) = 3.55*1e6; % NOK/km >300kV Duplex steel powerline cracle Veileder p142
% 
% DisCost = PlineC(PowerlineType(1))*D;

%% Cost Ontario case study
PlineC(1) = 198000; % $/km 25kV wood pole
PlineC(2) = 214000; % $/km 44kV wood pole
PlineC(3) = 553000; % $/km 115V steel frame
PlineC(4) = 1473000; % $/km 230kV steel frame
PlineC(5) = 1843000; % $/km 230kV steel frame double circuit

DisCost = PlineC(PowerlineType(1))*D;

if multifactor==1
    DisCost = DisCost * PWx;
end

end