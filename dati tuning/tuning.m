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
Q0 = Q(firstdate,:);
R0 = R(firstdate,:);
D0 = D(firstdate,:);
S0 = popolazione-E0-I0-Q(firstdate,:)-R(firstdate,:)-D(firstdate,:)-P0;

alpha=repmat(0.1:0.1:1,6,1)';
omega=repmat(0.1:0.1:1,6,1)';
for i =  1:size(alpha,1)
    for j = 1:size(omega,1)
        x0(7:12)=alpha(i,:);
        x0(13:18)=omega(j,:);
        [SEIQRDP_onNetwork{indexPeriod}.parametriR,SEIQRDP_onNetwork{indexPeriod}.resnormR,SEIQRDP_onNetwork{indexPeriod}.residualR] = fittingSEIQRDP_onNetwork(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,P0,E0,I0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
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
Q0 = Qmedia(firstdate,:);
R0 = Rmedia(firstdate,:);
D0 = Dmedia(firstdate,:);
S0 = popolazione-E0-I0-Qmedia(firstdate,:)-Rmedia(firstdate,:)-Dmedia(firstdate,:)-P0;

alpha=repmat(0.1:0.1:1,6,1)';
omega=repmat(0.1:0.1:1,6,1)';
for i =  1:size(alpha,1)
    for j = 1:size(omega,1)
        x0(7:12)=alpha(i,:);
        x0(13:18)=omega(j,:);
        [SEIQRDP_onNetwork{indexPeriod}.parametriR,SEIQRDP_onNetwork{indexPeriod}.resnormR,SEIQRDP_onNetwork{indexPeriod}.residualR] = fittingSEIQRDP_onNetwork(Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),popolazione,S0,P0,E0,I0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
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
Q0 = Qmedia(firstdate,:);
R0 = Rmedia(firstdate,:);
D0 = Dmedia(firstdate,:);
S0 = popolazione-E0-I0-Qmedia(firstdate,:)-Rmedia(firstdate,:)-Dmedia(firstdate,:)-P0;

alpha=repmat(0.1:0.1:1,6,1)';
omega=repmat(0.1:0.1:1,6,1)';
for i =  1:size(alpha,1)
    for j = 1:size(omega,1)
        x0(7:12)=alpha(i,:);
        x0(13:18)=omega(j,:);
        [SIQRDP_onNetwork{indexPeriod}.parametriR,SIQRDP_onNetwork{indexPeriod}.resnormR,SIQRDP_onNetwork{indexPeriod}.residualR] = fittingSIQRDP_onNetwork(Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),popolazione,S0,P0,I0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
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
Q0 = Q(firstdate,:);
R0 = R(firstdate,:);
D0 = D(firstdate,:);
S0 = popolazione-E0-I0-Q(firstdate,:)-R(firstdate,:)-D(firstdate,:)-P0;

alpha=repmat(0.1:0.1:1,6,1)';
omega=repmat(0.1:0.1:1,6,1)';
for i =  1:size(alpha,1)
    for j = 1:size(omega,1)
        x0(7:12)=alpha(i,:);
        x0(13:18)=omega(j,:);
        [SIQRDP_onNetwork{indexPeriod}.parametriR,SIQRDP_onNetwork{indexPeriod}.resnormR,SIQRDP_onNetwork{indexPeriod}.residualR] = fittingSIQRDP_onNetwork(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,P0,I0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
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

%save alpha.mat alpha
%save omega.mat omega
%save normErr_SEIQRDP_reale.mat normErr_SEIQRDP_reale
%save normErr_SEIQRDP_media.mat normErr_SEIQRDP_media
%save normErr_SIQRDP_reale.mat normErr_SIQRDP_reale
%save normErr_SIQRDP_media.mat normErr_SIQRDP_media

%%
time = unique(datetime(datestr(datenum(dataset{1}.data,'yyyy-mm-DDThh:MM:ss'))));
% If the number of confirmed Confirmed cases is small, it is difficult to know whether
% the quarantine has been rigorously applied or not. In addition, this
% suggests that the number of infectious is much larger than the number of
% confirmed cases
tStart = datetime(2020,03,22);
tEnd = datetime(2020,04,26);
indT = find(time>=tStart & time <=tEnd);

dataItalia = readtable("dati-regioni-xlsx/italia.xlsx");
timeItalia = unique(datetime(datestr(datenum(dataItalia.data,'yyyy-mm-DDThh:MM:ss'))));
load alpha.mat
load omega.mat
x0=[1.0*ones(1,6)]; %beta
x0=[x0  alpha(1,:)]; %alpha valori reali
x0=[x0  omega(1,:)]; %omega valori reali
x0=[x0 1/5*ones(1,6)]; %lambda -> periodo medio di incubazione è 5 giorni.
x0=[x0 0.5*ones(1,6)]; %delta
x0=[x0 1/11*ones(1,6)]; %gammaR

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
Q0 = Q(firstdate,:);
R0 = R(firstdate,:);
D0 = D(firstdate,:);
S0 = popolazione-E0-I0-Q(firstdate,:)-R(firstdate,:)-D(firstdate,:)-P0;

[SEIQRDP_onNetwork{indexPeriod}.parametriR,SEIQRDP_onNetwork{indexPeriod}.resnormR,SEIQRDP_onNetwork{indexPeriod}.residualR] = fittingSEIQRDP_onNetwork(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,P0,E0,I0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
SEIQRDP{indexPeriod}.YR = simulatedSEIQRDP_onNetwork(popolazione,S0,P0,E0,I0,Q0,R0,D0,Q(indT,:),R(indT,:),D(indT,:),L,SEIQRDP_onNetwork{indexPeriod}.parametriR,N,t,flag);
%plotSimulated(timeA,squeeze(SEIQRDP{indexPeriod}.YR(5,:,:)),squeeze(SEIQRDP{indexPeriod}.YR(6,:,:)),squeeze(SEIQRDP{indexPeriod}.YR(7,:,:)),time(indT),Tot(indT,:),Q(indT,:),R(indT,:),D(indT,:),'Dati reali: SEIQRDP on Network -',regioni)
x0(7:18)=[alpha(3,:) omega(3,:)];
[SIQRDP_onNetwork{indexPeriod}.parametriR,SIQRDP_onNetwork{indexPeriod}.resnormR,SIQRDP_onNetwork{indexPeriod}.residualR] = fittingSIQRDP_onNetwork(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,P0,I0,Q0,R0,D0,L,time(indT),x0([1:18 25:end]),flag,'Display','off','dt',0.1);
SIQRDP{indexPeriod}.YR= simulatedSIQRDP_onNetwork(popolazione,S0,P0,I0,Q0,R0,D0,Q(indT,:),R(indT,:),D(indT,:),L,SIQRDP_onNetwork{indexPeriod}.parametriR,N,t,flag);
%plotSimulated(timeA,squeeze(SIQRDP{indexPeriod}.YR(4,:,:)),squeeze(SIQRDP{indexPeriod}.YR(5,:,:)),squeeze(SIQRDP{indexPeriod}.YR(6,:,:)),time(indT),Tot(indT,:),Q(indT,:),R(indT,:),D(indT,:),'Dati reali: SIQRDP on Network -',regioni)

[SQRD_onNetwork{indexPeriod}.parametriR,SQRD_onNetwork{indexPeriod}.resnormR,SQRD_onNetwork{indexPeriod}.residualR] = fittingSQRD_onNetwork(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,Q0,R0,D0,L,time(indT),x0([1:6 31:end]),flag,'Display','off','dt',0.1);
SQRD{indexPeriod}.YR= simulatedSQRD_onNetwork(popolazione,S0,Q0,R0,D0,Q(indT,:),R(indT,:),D(indT,:),L,SQRD_onNetwork{indexPeriod}.parametriR,N,t,flag);
%plotSimulated(timeA,squeeze(SQRD{indexPeriod}.YR(2,:,:)),squeeze(SQRD{indexPeriod}.YR(3,:,:)),squeeze(SQRD{indexPeriod}.YR(4,:,:)),time(indT),Tot(indT,:),Q(indT,:),R(indT,:),D(indT,:),'Dati reali: SQRD on Network -',regioni)


x0(7:18)=[alpha(2,:) omega(2,:)];
x0(37:42)=rateDmedia(1,:);
I0 = 0.2*Qmedia(firstdate,:); % Initial number of infectious cases. Unknown but unlikely to be zero.
E0 = Qmedia(firstdate,:); % Initial number of exposed cases. Unknown but unlikely to be zero.
P0 = zeros(1,6);
Q0 = Qmedia(firstdate,:);
R0 = Rmedia(firstdate,:);
D0 = Dmedia(firstdate,:);
S0 = popolazione-E0-I0-Qmedia(firstdate,:)-Rmedia(firstdate,:)-Dmedia(firstdate,:)-P0;
[SEIQRDP_onNetwork{indexPeriod}.parametriM,SEIQRDP_onNetwork{indexPeriod}.resnormM,SEIQRDP_onNetwork{indexPeriod}.residualM] = fittingSEIQRDP_onNetwork(Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),popolazione,S0,P0,E0,I0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
SEIQRDP{indexPeriod}.YM = simulatedSEIQRDP_onNetwork(popolazione,S0,P0,E0,I0,Q0,R0,D0,Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),L,SEIQRDP_onNetwork{indexPeriod}.parametriM,N,t,flag);
%plotSimulated(timeA,squeeze(SEIQRDP{indexPeriod}.YM(5,:,:)),squeeze(SEIQRDP{indexPeriod}.YM(6,:,:)),squeeze(SEIQRDP{indexPeriod}.YM(7,:,:)),time(indT),Totmedia(indT,:),Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),'Dati in media mobile: SEIQRDP on Network -',regioni)

x0(7:18)=[alpha(4,:) omega(4,:)];
[SIQRDP_onNetwork{indexPeriod}.parametriM,SIQRDP_onNetwork{indexPeriod}.resnormM,SIQRDP_onNetwork{indexPeriod}.residualM] = fittingSIQRDP_onNetwork(Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),popolazione,S0,P0,I0,Q0,R0,D0,L,time(indT),x0([1:18 25:end]),flag,'Display','off','dt',0.1);
SIQRDP{indexPeriod}.YM= simulatedSIQRDP_onNetwork(popolazione,S0,P0,I0,Q0,R0,D0,Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),L,SIQRDP_onNetwork{indexPeriod}.parametriM,N,t,flag);
%plotSimulated(timeA,squeeze(SIQRDP{indexPeriod}.YM(4,:,:)),squeeze(SIQRDP{indexPeriod}.YM(5,:,:)),squeeze(SIQRDP{indexPeriod}.YM(6,:,:)),time(indT),Totmedia(indT,:),Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),'Dati in media mobile: SIQRDP on Network -',regioni)

[SQRD_onNetwork{indexPeriod}.parametriM,SQRD_onNetwork{indexPeriod}.resnormM,SQRD_onNetwork{indexPeriod}.residualM] = fittingSQRD_onNetwork(Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),popolazione,S0,Q0,R0,D0,L,time(indT),x0([1:6 31:end]),flag,'Display','off','dt',0.1);
SQRD{indexPeriod}.YM= simulatedSQRD_onNetwork(popolazione,S0,Q0,R0,D0,Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),L,SQRD_onNetwork{indexPeriod}.parametriM,N,t,flag);
%plotSimulated(timeA,squeeze(SQRD{indexPeriod}.YM(2,:,:)),squeeze(SQRD{indexPeriod}.YM(3,:,:)),squeeze(SQRD{indexPeriod}.YM(4,:,:)),time(indT),Totmedia(indT,:),Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),'Dati in media mobile: SQRD on Network -',regioni)

%save perido1.mat

%%

tStart = datetime(2020,04,27);
tEnd = datetime(2020,06,11);
indT = find(time>=tStart & time <=tEnd);
dt = 1/24; % time step
temp=time(indT);
timeA=datetime(temp(1), 'Locale', 'en_US'):dt:datetime(temp(end), 'Locale', 'en_US');
N = numel(timeA);
t = [0:N-1].*dt;

L=0.35*46*[0 117e03/popolazione(2) 0 108e03/popolazione(4) 0 0;
    117e03/popolazione(1) 0 302e03/popolazione(3) 300e03/popolazione(4) 0 0;
    0 302e03/popolazione(2) 0 136e03/popolazione(4) 0 0;
    108e03/popolazione(1) 300e03/popolazione(2) 136e03/popolazione(3) 0 311e03/popolazione(5)  116e03/popolazione(6);
    0 0 0 311e03/popolazione(4) 0 0;
    0 0 0 116e03/popolazione(4) 0 0];
Ldiag=diag(sum(L,1));
L=L-Ldiag;
indexPeriod=2;
flag=1; % diffusion yes
% Initial conditions

I0 = squeeze(SEIQRDP{indexPeriod-1}.YR(4,:,end)); 
E0 = squeeze(SEIQRDP{indexPeriod-1}.YR(3,:,end)); 
P0 = squeeze(SEIQRDP{indexPeriod-1}.YR(2,:,end));
S0 = squeeze(SEIQRDP{indexPeriod-1}.YR(1,:,end));
Q0 = squeeze(SEIQRDP{indexPeriod-1}.YR(5,:,end));
R0 = squeeze(SEIQRDP{indexPeriod-1}.YR(6,:,end));
D0 = squeeze(SEIQRDP{indexPeriod-1}.YR(7,:,end));

firstdate=find(time>=tStart & time<=tStart+datenum(1));
epsilonS= [0.0001 0.001 0.01 0.1 1];
epsilonI= [0.0001 0.001 0.01 0.1 1];
for i =  1:length(epsilonS)
    for j = 1:length(epsilonI)
        x0=[SEIQRDP_onNetwork{indexPeriod-1}.parametriR epsilonS(i) epsilonI(j)]; 
        [SEIQRDP_onNetwork{indexPeriod}.parametriR,SEIQRDP_onNetwork{indexPeriod}.resnormR,SEIQRDP_onNetwork{indexPeriod}.residualR] = fittingSEIQRDP_onNetwork(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,P0,E0,I0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
         for indR=1:6
           errR(1,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualR(1:length(indT),indR)));
           errR(2,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualR(length(indT)+2:2*length(indT),indR)));
           errR(3,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualR(2*length(indT)+2:3*length(indT),indR)));
         end
        normErr2_SEIQRDP_reale(i,j,:)=vecnorm(errR);
    end
end

for indR=1:6
    [x,indxEpsilonS] = min(normErr2_SEIQRDP_reale(:,:,indR));
    [y,indxEpsilonI] = min(x);
    epsilonS_SEIQRDP_reale(indR) = epsilonS(indxEpsilonS(indxEpsilonI));
    epsilonI_SEIQRDP_reale(indR) = epsilonI(indxEpsilonI);
end

%%
I0 = squeeze(SIQRDP{indexPeriod-1}.YR(3,:,end)); % Initial number of infectious cases. Unknown but unlikely to be zero.
P0 = squeeze(SIQRDP{indexPeriod-1}.YR(2,:,end));
S0 = squeeze(SIQRDP{indexPeriod-1}.YR(1,:,end));
Q0 = squeeze(SIQRDP{indexPeriod-1}.YR(4,:,end));
R0 = squeeze(SIQRDP{indexPeriod-1}.YR(5,:,end));
D0 = squeeze(SIQRDP{indexPeriod-1}.YR(6,:,end));

for i =  1:length(epsilonS)
    for j = 1:length(epsilonI)
        x0=[SIQRDP_onNetwork{indexPeriod-1}.parametriR epsilonS(i) epsilonI(j)]; 
        [SIQRDP_onNetwork{indexPeriod}.parametriR,SIQRDP_onNetwork{indexPeriod}.resnormR,SIQRDP_onNetwork{indexPeriod}.residualR] = fittingSIQRDP_onNetwork(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,P0,I0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
         for indR=1:6
           errR(1,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualR(1:length(indT),indR)));
           errR(2,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualR(length(indT)+2:2*length(indT),indR)));
           errR(3,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualR(2*length(indT)+2:3*length(indT),indR)));
         end
        normErr2_SIQRDP_reale(i,j,:)=vecnorm(errR);
    end
end

for indR=1:6
    [x,indxEpsilonS] = min(normErr2_SIQRDP_reale(:,:,indR));
    [y,indxEpsilonI] = min(x);
    epsilonS_SIQRDP_reale(indR) = epsilonS(indxEpsilonS(indxEpsilonI));
    epsilonI_SIQRDP_reale(indR) = epsilonI(indxEpsilonI);
end
%%

S0 = squeeze(SQRD{indexPeriod-1}.YR(1,:,end));
Q0 = squeeze(SQRD{indexPeriod-1}.YR(2,:,end));
R0 = squeeze(SQRD{indexPeriod-1}.YR(3,:,end));
D0 = squeeze(SQRD{indexPeriod-1}.YR(4,:,end));

epsilonS = [0.0001 0.001 0.01 0.1 1];
for i =  1:length(epsilonS)
     x0=[SQRD_onNetwork{indexPeriod-1}.parametriR epsilonS(i)];
     [SQRD_onNetwork{indexPeriod}.parametriR,SQRD_onNetwork{indexPeriod}.resnormR,SQRD_onNetwork{indexPeriod}.residualR] = fittingSQRD_onNetwork(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
     for indR=1:6
       errR(1,indR)=mean(abs(SQRD_onNetwork{indexPeriod}.residualR(1:length(indT),indR)));
       errR(2,indR)=mean(abs(SQRD_onNetwork{indexPeriod}.residualR(length(indT)+2:2*length(indT),indR)));
       errR(3,indR)=mean(abs(SQRD_onNetwork{indexPeriod}.residualR(2*length(indT)+2:3*length(indT),indR)));
     end
    normErr2_SQRD_reale(i,:)=vecnorm(errR);
end

for indR=1:6
    [x,indxEpsilonS] = min(normErr2_SQRD_reale(:,indR));
    epsilonS_SQRD_reale(indR) = epsilonS(indxEpsilonS);
end

%%
S0 = squeeze(SEIQRDP{indexPeriod-1}.YM(1,:,end));
P0 = squeeze(SEIQRDP{indexPeriod-1}.YM(2,:,end));
E0 = squeeze(SEIQRDP{indexPeriod-1}.YM(3,:,end));
I0 = squeeze(SEIQRDP{indexPeriod-1}.YM(4,:,end));
Q0 = squeeze(SEIQRDP{indexPeriod-1}.YM(5,:,end));
R0 = squeeze(SEIQRDP{indexPeriod-1}.YM(6,:,end));
D0 = squeeze(SEIQRDP{indexPeriod-1}.YM(7,:,end));

for i =  1:length(epsilonS)
    for j = 1:length(epsilonI)
        x0=[SEIQRDP_onNetwork{indexPeriod-1}.parametriM epsilonS(i) epsilonI(j)]; %% DA MODIFICARE
        [SEIQRDP_onNetwork{indexPeriod}.parametriM,SEIQRDP_onNetwork{indexPeriod}.resnormM,SEIQRDP_onNetwork{indexPeriod}.residualM] = fittingSEIQRDP_onNetwork(Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),popolazione,S0,P0,E0,I0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
         for indR=1:6
           errR(1,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualM(1:length(indT),indR)));
           errR(2,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualM(length(indT)+2:2*length(indT),indR)));
           errR(3,indR)=mean(abs(SEIQRDP_onNetwork{indexPeriod}.residualM(2*length(indT)+2:3*length(indT),indR)));
         end
        normErr2_SEIQRDP_media(i,j,:)=vecnorm(errR);
    end
end

for indR=1:6
    [x,indxEpsilonS] = min(normErr2_SEIQRDP_media(:,:,indR));
    [y,indxEpsilonI] = min(x);
    epsilonS_SEIQRDP_media(indR) = epsilonS(indxEpsilonS(indxEpsilonI));
    epsilonI_SEIQRDP_media(indR) = epsilonI(indxEpsilonI);
end

%%


S0 = squeeze(SIQRDP{indexPeriod-1}.YM(1,:,end));
P0 = squeeze(SIQRDP{indexPeriod-1}.YM(2,:,end));
I0 = squeeze(SIQRDP{indexPeriod-1}.YM(3,:,end));
Q0 = squeeze(SIQRDP{indexPeriod-1}.YM(4,:,end));
R0 = squeeze(SIQRDP{indexPeriod-1}.YM(5,:,end));
D0 = squeeze(SIQRDP{indexPeriod-1}.YM(6,:,end));
for i =  1:length(epsilonS)
    for j = 1:length(epsilonI)
        x0=[SIQRDP_onNetwork{indexPeriod-1}.parametriM epsilonS(i) epsilonI(j)]; %% DA MODIFICARE
        [SIQRDP_onNetwork{indexPeriod}.parametriM,SIQRDP_onNetwork{indexPeriod}.resnormM,SIQRDP_onNetwork{indexPeriod}.residualM] = fittingSIQRDP_onNetwork(Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),popolazione,S0,P0,I0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
         for indR=1:6
           errR(1,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualM(1:length(indT),indR)));
           errR(2,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualM(length(indT)+2:2*length(indT),indR)));
           errR(3,indR)=mean(abs(SIQRDP_onNetwork{indexPeriod}.residualM(2*length(indT)+2:3*length(indT),indR)));
         end
        normErr2_SIQRDP_media(i,j,:)=vecnorm(errR);
    end
end
for indR=1:6
    [x,indxEpsilonS] = min(normErr2_SIQRDP_media(:,:,indR));
    [y,indxEpsilonI] = min(x);
    epsilonS_SIQRDP_media(indR) = epsilonS(indxEpsilonS(indxEpsilonI));
    epsilonI_SIQRDP_media(indR) = epsilonI(indxEpsilonI);
end
%%
S0 = squeeze(SQRD{indexPeriod-1}.YM(1,:,end));
Q0 = squeeze(SQRD{indexPeriod-1}.YM(2,:,end));
R0 = squeeze(SQRD{indexPeriod-1}.YM(3,:,end));
D0 = squeeze(SQRD{indexPeriod-1}.YM(4 ,:,end));

epsilonS = [0.0001 0.001 0.01 0.1 1];
for i =  1:length(epsilonS)
    x0=[SQRD_onNetwork{indexPeriod-1}.parametriM epsilonS(i)];
    [SQRD_onNetwork{indexPeriod}.parametriM,SQRD_onNetwork{indexPeriod}.resnormM,SQRD_onNetwork{indexPeriod}.residualM] = fittingSQRD_onNetwork(Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),popolazione,S0,Q0,R0,D0,L,time(indT),x0,flag,'Display','off','dt',0.1);
     for indR=1:6
       errR(1,indR)=mean(abs(SQRD_onNetwork{indexPeriod}.residualM(1:length(indT),indR)));
       errR(2,indR)=mean(abs(SQRD_onNetwork{indexPeriod}.residualM(length(indT)+2:2*length(indT),indR)));
       errR(3,indR)=mean(abs(SQRD_onNetwork{indexPeriod}.residualM(2*length(indT)+2:3*length(indT),indR)));
     end
    normErr2_SQRD_media(i,:)=vecnorm(errR);
end

for indR=1:6
    [x,indxEpsilonS] = min(normErr2_SQRD_media(:,indR));
    epsilonS_SQRD_media(indR) = epsilonS(indxEpsilonS);
end

%%

epsilonS = [epsilonS_SEIQRDP_reale; epsilonS_SEIQRDP_media ; epsilonS_SIQRDP_reale ;epsilonS_SIQRDP_media; epsilonS_SQRD_reale ;epsilonS_SQRD_media];
%save epsilonS.mat epsilonS

epsilonI = [epsilonI_SEIQRDP_reale; epsilonI_SEIQRDP_media ; epsilonI_SIQRDP_reale ;epsilonI_SIQRDP_media];
%save epsilonI.mat epsilonI