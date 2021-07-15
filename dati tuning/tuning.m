%% pre-processing

clc
clear 
close all

numRegioni=6;
regioni={'Piemonte','Lombardia','Veneto','Emilia-Romangna','Marche','Toscana'};

dataset{1} = readtable("dati-regioni-xlsx/piemonte.xlsx");
dataset{2} = readtable("dati-regioni-xlsx/lombardia.xlsx");
dataset{3} = readtable("dati-regioni-xlsx/veneto.xlsx");
dataset{4} = readtable("dati-regioni-xlsx/emilia.xlsx");
dataset{5} = readtable("dati-regioni-xlsx/marche.xlsx");
dataset{6} = readtable("dati-regioni-xlsx/toscana.xlsx");

time = unique(datetime(datestr(datenum(dataset{1}.data,'yyyy-mm-DDThh:MM:ss'))));

tStart = datetime(2020,03,22);
periodoAnalisi = find(time>=tStart);

dataset{1}.dimessi_guariti(152)=26669;
dataset{1}.totale_casi(296)=186172;
dataset{1}.totale_casi(297)=186387;
dataset{1}.totale_casi(298)=187270;
dataset{3}.dimessi_guariti(119)=16670;
dataset{3}.dimessi_guariti(120)=16679;
dataset{3}.dimessi_guariti(155)=17025;
dataset{3}.deceduti(239)=2265;
dataset{5}.dimessi_guariti(295)=21756;
dataset{5}.deceduti(60:95)=dataset{5}.deceduti(60:95)-11;
dataset{5}.deceduti(90:122)=dataset{5}.deceduti(90:122)-3;
dataset{5}.deceduti(73:128)=dataset{5}.deceduti(73:128)-4;
dataset{5}.totale_casi(85)=6668;
dataset{5}.totale_casi(146)=6810;


R=zeros(length(dataset{1}.data),numRegioni);
Tot=zeros(length(dataset{1}.data),numRegioni);
D=zeros(length(dataset{1}.data),numRegioni);
Tamp=zeros(length(dataset{1}.data),numRegioni);

tStart = datetime(2020,03,22);
periodoAnalisi = find(time>=tStart);

Rmedia=zeros(length(dataset{1}.data),numRegioni);
Dmedia=zeros(length(dataset{1}.data),numRegioni);
Totmedia=zeros(length(dataset{1}.data),numRegioni);
Tampmedia=zeros(length(dataset{1}.data),numRegioni);
for indR=1:numRegioni
    Tot(:,indR)=dataset{indR}.totale_casi;
    R(:,indR)=dataset{indR}.dimessi_guariti;
    D(:,indR)=dataset{indR}.deceduti;
    Tamp(:,indR)=dataset{indR}.tamponi;
    Rmedia(:,indR)=smooth(R(:,indR),7); %moving average
    Totmedia(:,indR)=smooth(Tot(:,indR),7);
    Dmedia(:,indR)=smooth(D(:,indR),7);
end
Q=Tot-R-D;
Qmedia=Totmedia-Rmedia-Dmedia;

%% tuning alpha e omega
x0=[1.0*ones(1,6)]; %beta
x0=[x0  0.4 0.7 0.4 0.8 0.5 0.5]; %alpha valori reali
x0=[x0  0.1 0.5 0.3 0.2 0.5 0.5]; %omega valori reali
x0=[x0 1/5*ones(1,6)]; %lambda -> periodo medio di incubazione è 5 giorni.
x0=[x0 0.5*ones(1,6)]; %delta
x0=[x0 1/11*ones(1,6)]; %gammaR

time = unique(datetime(datestr(datenum(dataset{1}.data,'yyyy-mm-DDThh:MM:ss'))));
% If the number of confirmed Confirmed cases is small, it is difficult to know whether
% the quarantine has been rigorously applied or not. In addition, this
% suggests that the number of infectious is much larger than the number of
% confirmed cases
tStart = datetime(2020,03,22);
tEnd = datetime(2020,04,26);
indT = find(time>=tStart & time <=tEnd);

% tasso di morte
rateD=[];
rateDmedia=[];
for indR=1:6
    rateD(:,indR) = (diff(D(:,indR))./diff(datenum(time-time(1))))./Q(2:end,indR);
    rateDmedia(:,indR) = (diff(Dmedia(:,indR))./diff(datenum(time-time(1))))./Qmedia(2:end,indR);
end
rateD=rateD(periodoAnalisi-1,:);
rateDmedia=rateDmedia(periodoAnalisi-1,:);

popolazione=[4.311e6 10.028e6 4.879e6 4.464e6 1.513e6 3.693e6]; 

L=zeros(6,6);
indexPeriod=1;
flag=0;
% Initial conditions
dt = 1/24; % time step
temp=time(indT);
timeA=datetime(temp(1), 'Locale', 'en_US'):dt:datetime(temp(end), 'Locale', 'en_US');
N = numel(timeA);
t = [0:N-1].*dt;

firstdate=find(time>=tStart & time<=tStart+datenum(1));
x0 = [x0 rateD(1,:)]; %gammaD
I0 = 0.2*Q(firstdate,:); % Initial number of infectious cases. Unknown but unlikely to be zero.
E0 = Q(firstdate,:); % Initial number of exposed cases. Unknown but unlikely to be zero.
P0 = zeros(1,6);
S0 = popolazione-E0-I0-Q(firstdate,:)-R(firstdate,:)-D(firstdate,:)-P0;

alpha=repmat(0.1:0.1:1,6,1)';
omega=repmat(0.1:0.1:1,6,1)';
for i =  1:size(alpha,1)
    for j = 1:size(omega,1)
        i
        j
        x0(7:12)=alpha(i,:);
        x0(13:18)=omega(j,:);
        [SEIQRDP_onNetwork{indexPeriod}.parametriR,SEIQRDP_onNetwork{indexPeriod}.resnormR,SEIQRDP_onNetwork{indexPeriod}.residualR] = fittingSEIQRDP_onNetwork(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,E0,I0,P0,L,time(indT),x0,flag,'Display','off','dt',0.1);
         for indR=1:6
           errR(1,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualR(1:length(indT),indR)));
           errR(2,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualR(length(indT)+2:2*length(indT),indR)));
           errR(3,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualR(2*length(indT)+2:3*length(indT),indR)));
         end
        normErrR(i,j,:)=vecnorm(errR);
    end
end
%%
for indR=1:6
    [x,indxAlpha] = min(normErrR(:,:,indR));
    [y,indxOmega] = min(x);
    alphaSEQIRDP_reale(indR) = alpha(indxAlpha(indxOmega),1);
    omegaSEQIRDP_reale(indR) = omega(indxOmega,1); 
end

%save SEIQRDP_reale.mat

%%
x0=[1.0*ones(1,6)]; %beta
x0=[x0  0.4 0.7 0.4 0.8 0.5 0.5]; %alpha valori reali
x0=[x0  0.1 0.5 0.3 0.2 0.5 0.5]; %omega valori reali
x0=[x0 1/5*ones(1,6)]; %lambda -> periodo medio di incubazione è 5 giorni.
x0=[x0 0.5*ones(1,6)]; %delta
x0=[x0 1/11*ones(1,6)]; %gammaR
firstdate=find(time>=tStart & time<=tStart+datenum(1));
x0 = [x0 rateDmedia(1,:)]; %gammaD
I0 = 0.2*Qmedia(firstdate,:); % Initial number of infectious cases. Unknown but unlikely to be zero.
E0 = Qmedia(firstdate,:); % Initial number of exposed cases. Unknown but unlikely to be zero.
P0 = zeros(1,6);
S0 = popolazione-E0-I0-Qmedia(firstdate,:)-Rmedia(firstdate,:)-Dmedia(firstdate,:)-P0;

alpha=repmat(0.1:0.1:1,6,1)';
omega=repmat(0.1:0.1:1,6,1)';
for i =  1:size(alpha,1)
    for j = 1:size(omega,1)
        i
        j
        x0(7:12)=alpha(i,:);
        x0(13:18)=omega(j,:);
        [SEIQRDP_onNetwork{indexPeriod}.parametriR,SEIQRDP_onNetwork{indexPeriod}.resnormR,SEIQRDP_onNetwork{indexPeriod}.residualR] = fittingSEIQRDP_onNetwork(Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),popolazione,S0,E0,I0,P0,L,time(indT),x0,flag,'Display','off','dt',0.1);
         for indR=1:6
           errR(1,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualR(1:length(indT),indR)));
           errR(2,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualR(length(indT)+2:2*length(indT),indR)));
           errR(3,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualR(2*length(indT)+2:3*length(indT),indR)));
         end
        normErr_SEIQRDP_media(i,j,:)=vecnorm(errR);
    end
end
%%
for indR=1:6
    [x,indxAlpha] = min(normErr_SEIQRDP_media(:,:,indR));
    [y,indxOmega] = min(x);
    alphaSEIQRDP_media(indR) = alpha(indxAlpha(indxOmega),1);
    omegaSEIQRDP_media(indR) = omega(indxOmega,1);
end

%save SEIQRDP_media.mat

%%
x0=[1.0*ones(1,6)]; %beta
x0=[x0  0.4 0.7 0.4 0.8 0.5 0.5]; %alpha valori reali
x0=[x0  0.1 0.5 0.3 0.2 0.5 0.5]; %omega valori reali
x0=[x0 0.5*ones(1,6)]; %delta
x0=[x0 1/11*ones(1,6)]; %gammaR
firstdate=find(time>=tStart & time<=tStart+datenum(1));
x0 = [x0 rateDmedia(1,:)]; %gammaD
I0 = 0.2*Qmedia(firstdate,:); % Initial number of infectious cases. Unknown but unlikely to be zero.
P0 = zeros(1,6);
S0 = popolazione-E0-I0-Qmedia(firstdate,:)-Rmedia(firstdate,:)-Dmedia(firstdate,:)-P0;

alpha=repmat(0.1:0.1:1,6,1)';
omega=repmat(0.1:0.1:1,6,1)';
for i =  1:size(alpha,1)
    for j = 1:size(omega,1)
        i
        j
        x0(7:12)=alpha(i,:);
        x0(13:18)=omega(j,:);
        [SIQRDP_onNetwork{indexPeriod}.parametriR,SIQRDP_onNetwork{indexPeriod}.resnormR,SIQRDP_onNetwork{indexPeriod}.residualR] = fittingSIQRDP_onNetwork(Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),popolazione,S0,I0,P0,L,time(indT),x0,flag,'Display','off','dt',0.1);
         for indR=1:6
           errR(1,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualR(1:length(indT),indR)));
           errR(2,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualR(length(indT)+2:2*length(indT),indR)));
           errR(3,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualR(2*length(indT)+2:3*length(indT),indR)));
         end
        normErr_SIQRDP_media(i,j,:)=vecnorm(errR);
    end
end

for indR=1:6
    [x,indxAlpha] = min(normErr_SIQRDP_media(:,:,indR));
    [y,indxOmega] = min(x);
    alphaSIQRDP_media(indR) = alpha(indxAlpha(indxOmega),1);
    omegaSIQRDP_media(indR) = omega(indxOmega,1);
end

save SIQRDP_media.mat

%%
x0=[1.0*ones(1,6)]; %beta
x0=[x0  0.4 0.7 0.4 0.8 0.5 0.5]; %alpha valori reali
x0=[x0  0.1 0.5 0.3 0.2 0.5 0.5]; %omega valori reali
x0=[x0 0.5*ones(1,6)]; %delta
x0=[x0 1/11*ones(1,6)]; %gammaR
firstdate=find(time>=tStart & time<=tStart+datenum(1));
x0 = [x0 rateD(1,:)]; %gammaD
I0 = 0.2*Q(firstdate,:); % Initial number of infectious cases. Unknown but unlikely to be zero.
P0 = zeros(1,6);
S0 = popolazione-E0-I0-Q(firstdate,:)-R(firstdate,:)-D(firstdate,:)-P0;

alpha=repmat(0.1:0.1:1,6,1)';
omega=repmat(0.1:0.1:1,6,1)';
for i =  1:size(alpha,1)
    for j = 1:size(omega,1)
        i
        j
        x0(7:12)=alpha(i,:);
        x0(13:18)=omega(j,:);
        [SIQRDP_onNetwork{indexPeriod}.parametriR,SIQRDP_onNetwork{indexPeriod}.resnormR,SIQRDP_onNetwork{indexPeriod}.residualR] = fittingSIQRDP_onNetwork(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,I0,P0,L,time(indT),x0,flag,'Display','off','dt',0.1);
         for indR=1:6
           errR(1,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualR(1:length(indT),indR)));
           errR(2,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualR(length(indT)+2:2*length(indT),indR)));
           errR(3,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualR(2*length(indT)+2:3*length(indT),indR)));
         end
        normErr_SIQRDP_reale(i,j,:)=vecnorm(errR);
    end
end

for indR=1:6
    [x,indxAlpha] = min(normErr_SIQRDP_reale(:,:,indR));
    [y,indxOmega] = min(x);
    alphaSIQRDP_reale(indR) = alpha(indxAlpha(indxOmega),1);
    omegaSIQRDP_reale(indR) = omega(indxOmega,1);
end

save SIQRDP_reale.mat

%%

alpha =[alphaSEIQRDP_reale; alphaSEIQRDP_media; alphaSIQRDP_reale; alphaSIQRDP_media];
omega = [omegaSEIQRDP_reale;omegaSEIQRDP_media; omegaSIQRDP_reale;omegaSIQRDP_media];

save alpha.mat alpha
save omega.mat omega
save normErr_SEIQRDP_reale.mat normErr_SEIQRDP_reale
save normErr_SEIQRDP_media.mat normErr_SEIQRDP_media
save normErr_SIQRDP_reale.mat normErr_SIQRDP_reale
save normErr_SIQRDP_media.mat normErr_SIQRDP_media