function [Y] = RK4_onNetwork(Fun,Y,H,L,Z,dt)
    % NUmerical trick: the parameters are assumed constant between
    % two time steps
    % Runge-Kutta of order 4
    k_1 = Fun(Y,H,L,Z);
    k_2 = Fun(Y+0.5*dt*k_1,H,L,Z);
    k_3 = Fun(Y+0.5*dt*k_2,H,L,Z);
    k_4 = Fun(Y+k_3*dt,H,L,Z);
    % output
    Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
end


