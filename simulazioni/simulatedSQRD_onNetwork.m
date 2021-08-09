function [Y] = simulatedSQRD_onNetwork(Npop,S0,Q0,R0,D0,Q,R,D,L,para,N,t,flag)
        % I simply rename the inputs
        beta = abs(para(1:6)); %rimane invariato
        gammaR = abs(para(7:12));
        gammaD = abs(para(13:18));
        if flag==1
            epsilonS=abs(para(19));
        else
            epsilonS=0;
        end
        Y = zeros(4,6,N); %  There are seven different states for each region
        for j=1:6
            Y(1,j,1) = S0(j);
            Y(2,j,1) = Q0(j);
            Y(3,j,1) = R0(j);
            Y(4,j,1) = D0(j);
        end
    
    dt = median(diff(t));
    
    
        modelFun = @(Y,H,L,Z) H*Y*L' + Z;
        
        % Very large recovery rate should not occur but can lead to
        % numerical errors.
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
    
    % Y = round(Y);
end

