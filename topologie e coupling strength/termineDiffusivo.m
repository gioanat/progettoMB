%% influenza del termine diffusivo

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

tStart = datetime(2020,03,10);
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

rateD=[];
rateDmedia=[];
for indR=1:6
    rateD(:,indR) = (diff(D(:,indR))./diff(datenum(time-time(1))))./Q(2:end,indR);
    rateDmedia(:,indR) = (diff(Dmedia(:,indR))./diff(datenum(time-time(1))))./Qmedia(2:end,indR);
end
rateD=rateD(periodoAnalisi-1,:);
rateDmedia=rateDmedia(periodoAnalisi-1,:);

popolazione=[4.311e6 10.028e6 4.879e6 4.464e6 1.513e6 3.693e6]; 

time = unique(datetime(datestr(datenum(dataset{1}.data,'yyyy-mm-DDThh:MM:ss'))));
% If the number of confirmed Confirmed cases is small, it is difficult to know whether
% the quarantine has been rigorously applied or not. In addition, this
% suggests that the number of infectious is much larger than the number of
% confirmed cases
tEnd = datetime(2020,07,01);
indT = find(time>=tStart & time <=tEnd);
load alpha
load omega
x0=[1.0*ones(1,6)]; %beta
x0=[x0  alpha(1,:)]; %alpha valori reali
x0=[x0  omega(1,:)]; %omega valori reali
x0=[x0 0.5*ones(1,6)]; %delta
x0=[x0 1/5*ones(1,6)]; %lambda -> periodo medio di incubazione è 5 giorni.
x0=[x0 1/11*ones(1,6)]; %gammaR
firstdate=find(time>=tStart & time<=tStart+datenum(1));
x0 = [x0 rateD(1,:)]; %gammaD

indT = find(time>=tStart & time <=tEnd);
dt = 1/24; % time step
temp=time(indT);
timeA=datetime(temp(1), 'Locale', 'en_US'):dt:datetime(temp(end), 'Locale', 'en_US');
N = numel(timeA);
t = [0:N-1].*dt;

L=0.35*60*[0 117e03/popolazione(2) 0 108e03/popolazione(4) 0 0;
    117e03/popolazione(1) 0 302e03/popolazione(3) 300e03/popolazione(4) 0 0;
    0 302e03/popolazione(2) 0 136e03/popolazione(4) 0 0;
    108e03/popolazione(1) 300e03/popolazione(2) 136e03/popolazione(3) 0 311e03/popolazione(5)  116e03/popolazione(6);
    0 0 0 311e03/popolazione(4) 0 0;
    0 0 0 116e03/popolazione(4) 0 0];
% L=ones(6,6);
% L=L-diag(diag(L));
Ldiag=diag(sum(L,1));
L=L-Ldiag;
indexPeriod=1;
flag=1; % diffusion yes

I0 = 0.2*Q(firstdate,:); % Initial number of infectious cases. Unknown but unlikely to be zero.
E0 = Q(firstdate,:); % Initial number of exposed cases. Unknown but unlikely to be zero.
P0 = zeros(1,6);
Q0 = Q(firstdate,:);
R0 = R(firstdate,:);
D0 = D(firstdate,:);
S0 = popolazione-E0-I0-Q(firstdate,:)-R(firstdate,:)-D(firstdate,:)-P0;

timeFitting=[];
temp=time(indT);
timeFitting=[timeFitting datetime(temp(1), 'Locale', 'en_US'):dt:datetime(temp(end), 'Locale', 'en_US')];

epsilonI = [10^(-7) 1];
epsilonS = [10^(-7) 1];
for j = 1:length(epsilonI)
    [SEIQRDP_onNetwork{indexPeriod}.parametriR,SEIQRDP_onNetwork{indexPeriod}.resnormR,SEIQRDP_onNetwork{indexPeriod}.residualR] = fittingSEIQRDP_termDiff(Q(indT,:),R(indT,:),D(indT,:),popolazione,S0,P0,E0,I0,Q0,R0,D0,L,time(indT),x0,epsilonS(j),epsilonI(j),flag,'Display','off','dt',0.1);
    SEIQRDP{indexPeriod}.YR = simulatedSEIQRDP_onNetwork(popolazione,S0,P0,E0,I0,Q0,R0,D0,Q(indT,:),R(indT,:),D(indT,:),L,[SEIQRDP_onNetwork{indexPeriod}.parametriR epsilonS(j) epsilonI(j)],N,t,flag);
    figure
    hold on
    plot(timeFitting, squeeze(SEIQRDP{indexPeriod}.YR(5,1,:)),'k',timeFitting, squeeze(SEIQRDP{indexPeriod}.YR(5,2,:)),'m',timeFitting, squeeze(SEIQRDP{indexPeriod}.YR(5,3,:)),'c',timeFitting, squeeze(SEIQRDP{indexPeriod}.YR(5,4,:)),'r',timeFitting, squeeze(SEIQRDP{indexPeriod}.YR(5,5,:)),'g',timeFitting, squeeze(SEIQRDP{indexPeriod}.YR(5,6,:)),'b')
    y = ylim; % current y-axis limits
    [x,indx]=max(squeeze(SEIQRDP{indexPeriod}.YR(5,1,:)));
    plot([timeFitting(indx) timeFitting(indx)],[y(1) y(2)],'k:',timeFitting(indx),x,'k*')
    [x,indx]=max(squeeze(SEIQRDP{indexPeriod}.YR(5,2,:)));
    plot([timeFitting(indx) timeFitting(indx)],[y(1) y(2)],'m:',timeFitting(indx),x,'m*')
    [x,indx]=max(squeeze(SEIQRDP{indexPeriod}.YR(5,3,:)));
    plot([timeFitting(indx) timeFitting(indx)],[y(1) y(2)],'c:',timeFitting(indx),x,'c*')
    [x,indx]=max(squeeze(SEIQRDP{indexPeriod}.YR(5,4,:)));
    plot([timeFitting(indx) timeFitting(indx)],[y(1) y(2)],'r:',timeFitting(indx),x,'r*')
    [x,indx]=max(squeeze(SEIQRDP{indexPeriod}.YR(5,5,:)));
    plot([timeFitting(indx) timeFitting(indx)],[y(1) y(2)],'g:',timeFitting(indx),x,'g*')
    [x,indx]=max(squeeze(SEIQRDP{indexPeriod}.YR(5,6,:)));
    plot([timeFitting(indx) timeFitting(indx)],[y(1) y(2)],'b:',timeFitting(indx),x,'b*')
    legend('Piemonte','Lombardia','Veneto','Emilia-Romagna','Marche','Toscana')
    xlabel("time")
    ylabel("Q")
    if j==1
       title("\epsilon_I, \epsilon_S = 10^{-7}")
    else
       title("\epsilon_I, \epsilon_S = 1")
    end
    hold off
end

function [x,resnorm,residual] = fittingSEIQRDP_termDiff(Q,R,D,Npop,S0,P0,E0,I0,Q0,R0,D0,L,time,x0,epsilonS,epsilonI,flag,varargin)
% per stimare i parametri utilizziamo lsqcurvefit (stimatore dei minimi
% guadrati -> minimizza i residui) il quale prende in input
% - fun
% - x0: stima iniziale dei parametri
% - xdata: tempo
% - ydata: i valori noti
% - lb and ub: limite inferiore e limite superiore dei parametri

% - options
%Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('tolX',1e-05);  %  option for optimset
p.addOptional('tolFun',1e-05);  %  option for optimset
p.addOptional('Display','iter'); % Display option for optimset
p.addOptional('dt',0.1); % time step for the fitting

p.parse(varargin{:});
tolX = p.Results.tolX ;
tolFun = p.Results.tolFun ;
Display  = p.Results.Display ;
dt  = p.Results.dt ;

% Options for lsqcurvfit
options=optimset('TolX',tolX,'TolFun',tolFun,...
    'MaxFunEvals',1200,'Display',Display);
fs = 1./dt;
tTarget = round(datenum(time-time(1))*fs)/fs; % Number of days with one decimal
t = tTarget(1):dt:tTarget(end); % oversample to ensure that the algorithm converges
ydata = [Q;R;D];
if flag==1
    ub = ones(1,42); % upper bound of the parameters
    lb = zeros(1,42); % lower bound of the parameters
else
    ub = ones(1,42); % upper bound of the parameters
    lb = zeros(1,42); % lower bound of the parameters
end

[x,resnorm,residual]= lsqcurvefit(@(para,t) fun(para,t),x0,tTarget,ydata,lb,ub,options);

function [output] = fun(para,t0)

        % I simply rename the inputs
        beta = abs(para(1:6)); %rimane invariato
        alpha = abs(para(7:12));
        omega = abs(para(13:18));
        lambda = abs(para(19:24));
        delta = abs(para(25:30));
        gammaR = abs(para(31:36));
        gammaD = abs(para(37:42));
        
        % Initial conditions
        N = numel(t);
        Y = zeros(7,6,N); %  There are seven different states for each region
        for j=1:6
            Y(1,j,1) = S0(j);
            Y(2,j,1) = P0(j);
            Y(3,j,1) = E0(j);
            Y(4,j,1) = I0(j);
            Y(5,j,1) = Q0(j);
            Y(6,j,1) = R0(j);
            Y(7,j,1) = D0(j);
        end

        modelFun = @(Y,H,L,Z) H*Y*L' + Z;
        
        % Very large recovery rate should not occur but can lead to
        % numerical errors.
        %if lambda>10, warning('lambda is abnormally high'); end
        
        % ODE resolution
        for ii=1:N-1
            Z=zeros(7,6);
            for j=1:6
                A = [-alpha(j) omega(j) 0 0 0 0 0;
                    alpha(j) -omega(j) 0 0 0 0 0;
                    0 0 -lambda(j) 0 0 0 0;
                    0 0 lambda(j) -delta(j) 0 0 0;
                    0 0 0 delta(j) -(gammaR(j)+gammaD(j)) 0 0;
                    0 0 0 0 gammaR(j) 0 0;
                    0 0 0 0 gammaD(j) 0 0];
%                 A = [0 omega(j) 0 0 0 0 0;
%                     0 -omega(j) 0 0 0 0 0;
%                     0 0 -lambda(j) 0 0 0 0;
%                     0 0 lambda(j) -delta(j) 0 0 0;
%                     0 0 0 delta(j) -(gammaR(j)+gammaD(j)) 0 0;
%                     0 0 0 0 gammaR(j) 0 0;
%                     0 0 0 0 gammaD(j) 0 0];
                SI = Y(1,j,ii)*Y(4,j,ii);
                F = zeros(7,1);
                F([1 3],1) = [-beta(j)/Npop(j);beta(j)/Npop(j)].*SI;
%                 SI = Y(1,j,ii)*Y(5,j,ii);
%                 F([1 2],1) = [-alpha(j)/Npop(j);alpha(j)/Npop(j)].*SI; 
                Z(:,j)=A*Y(:,j,ii)+F;
            end
            H=diag([epsilonS 0 0 epsilonI 0 0 0]);
            Y(:,:,ii+1) = RK4_onNetwork(modelFun,Y(:,:,ii),H,L,Z,dt);
        end
        for j=1:6
            Q1(:,j) = Y(5,j,1:N);
            R1(:,j) = Y(6,j,1:N);
            D1(:,j) = Y(7,j,1:N);

            Q2(:,j) = interp1(t,Q1(:,j),t0);
            R2(:,j) = interp1(t,R1(:,j),t0);
            D2(:,j) = interp1(t,D1(:,j),t0);
        end
        output = ([Q2;R2;D2]);
end

end