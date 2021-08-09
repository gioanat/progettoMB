function [Y] = simulatedSIQRDP_onNetwork(Npop,S0,P0,I0,Q0,R0,D0,Q,R,D,L,para,N,t,flag)
        % I simply rename the inputs
        beta = abs(para(1:6)); %rimane invariato
        alpha = abs(para(7:12));
        omega = abs(para(13:18));
        delta = abs(para(19:24));
        gammaR = abs(para(25:30));
        gammaD = abs(para(31:36));
        if flag==1
            epsilonS=abs(para(37));
            epsilonI=abs(para(38));
        else
            epsilonS=0;
            epsilonI=0;
        end
        Y = zeros(6,6,N); %  There are seven different states for each region
        for j=1:6
            Y(1,j,1) = S0(j);
            Y(2,j,1) = P0(j);
            Y(3,j,1) = I0(j);
            Y(4,j,1) = Q0(j);
            Y(5,j,1) = R0(j);
            Y(6,j,1) = D0(j);
        end
    
    dt = median(diff(t));
    
    
        modelFun = @(Y,H,L,Z) H*Y*L' + Z;
        
        % Very large recovery rate should not occur but can lead to
        % numerical errors.
        for ii=1:N-1
            Z=zeros(6,6);
            for j=1:6
                Npop(j) = sum(squeeze(Y(:,j,ii)));
                A = [-alpha(j) omega(j) 0 0 0 0;
                    alpha(j) -omega(j) 0 0 0 0;
                    0 0 -delta(j) 0 0 0;
                    0 0 delta(j) -(gammaR(j)+gammaD(j)) 0 0;
                    0 0 0 gammaR(j) 0 0;
                    0 0 0 gammaD(j) 0 0];
                SI = Y(1,j,ii)*Y(3,j,ii);
                F = zeros(6,1);
                F([1 3],1) = [-beta(j)/Npop(j);beta(j)/Npop(j)].*SI;
                Z(:,j)=A*Y(:,j,ii)+F;
            end
            H=diag([epsilonS 0 epsilonI 0 0 0]);
            Y(:,:,ii+1) = RK4_onNetwork(modelFun,Y(:,:,ii),H,L,Z,dt);
        end
    
    % Y = round(Y);
end