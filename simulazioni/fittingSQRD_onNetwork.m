function [x,resnorm,residual] = fittingSQRD_onNetwork(Q,R,D,Npop,S0,Q0,R0,D0,L,time,x0,flag,varargin)
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
    ub = [1*ones(1,6) ones(1,13)]; % upper bound of the parameters
    lb = zeros(1,19); % lower bound of the parameters
else
    ub = [1*ones(1,6) ones(1,12)]; % upper bound of the parameters
    lb = zeros(1,18); % lower bound of the parameters
end

[x,resnorm,residual]= lsqcurvefit(@(para,t) fun(para,t),x0,tTarget,ydata,lb,ub,options);

function [output] = fun(para,t0)

        % I simply rename the inputs
        beta = abs(para(1:6));
        gammaR = abs(para(7:12));
        gammaD = abs(para(13:18));
        if flag==1
            epsilonS=abs(para(19));
        else
            epsilonS=0;
        end
        
        % Initial conditions
        N = numel(t);
        Y = zeros(4,6,N); %  There are seven different states for each region
        for j=1:6
            Y(1,j,1) = S0(j);
            Y(2,j,1) = Q0(j);
            Y(3,j,1) = R0(j);
            Y(4,j,1) = D0(j);
        end

        modelFun = @(Y,H,L,Z) H*Y*L' + Z;
        
        % Very large recovery rate should not occur but can lead to
        % numerical errors.
        %if lambda>10, warning('lambda is abnormally high'); end
        
        % ODE resolution
        for ii=1:N-1
            Z=zeros(4,6);
            for j=1:6
                Npop(j) = sum(squeeze(Y(:,j,ii)));
                A = [0 0 0 0;
                    0 -(gammaR(j)+gammaD(j)) 0 0;
                    0 gammaR(j) 0 0;
                    0 gammaD(j) 0 0];
                SI = Y(1,j,ii)*Y(2,j,ii);
                F = zeros(4,1);
                F([1 2],1) = [-beta(j)/Npop(j);beta(j)/Npop(j)].*SI;
                Z(:,j)=A*Y(:,j,ii)+F;
            end
            H=diag([epsilonS 0 0 0]);
            Y(:,:,ii+1) = RK4_onNetwork(modelFun,Y(:,:,ii),H,L,Z,dt);
        end
        for j=1:6
            Q1(:,j) = Y(2,j,1:N);
            R1(:,j) = Y(3,j,1:N);
            D1(:,j) = Y(4,j,1:N);

            Q2(:,j) = interp1(t,Q1(:,j),t0);
            R2(:,j) = interp1(t,R1(:,j),t0);
            D2(:,j) = interp1(t,D1(:,j),t0);
        end
        output = ([Q2;R2;D2]);
end

end