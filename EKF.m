%       title: EKF.m
%      author: Elijah Vautour
%        date: November 6, 2022
% description: The EKF class contains the functions used in the reduced
%              dynamic extended Kalman filter (RDEKF) algorithm and generates a
%              configurable EKF filter object which stores user-defined filter
%              properties.

classdef EKF < matlab.mixin.Copyable
    % calculated properties
    properties (Hidden)
        Q % covariance of process noise
    end
    
    % set properties
    properties
        % EKF timestep parameters
        dt      %   [s], simulation time step
        dtMes   %   [s], steady state measurement time step
        propSet % [min], propagation time in between measurement updates
        mesSet  % [min], duration of measurement update phase
        % EKF model definition parameters
        proMode % [-], process model mode
                %       ->   'TB' = simple two-body model
                %       -> 'TBJ2' = two-body model with J2 perturbation
        mesMode % [-], measurement model mode
                %       ->  'PR' = pseudorange observables
                %       -> 'SPS' = receiver computed eci observables
        % process model variance parameters
        proVar = struct('P',[], ...  %     [m^2], process model: position variance
                        'V',[], ...  % [m^2/s^2], process model: velocity variance
                        'Sf',[], ... %       [?], process model: clock variance parameter #1
                        'Sb',[]);    %       [?], process model: clock variance parameter #2
        % measurement model variance parameters
        mesVar = struct('Pr',[], ... %     [m^2], measurement model: pseudorange variance
                        'P',[], ...  %     [m^2], measurement model: position variance
                        'V',[]);     % [m^2/s^2], measurement model: velocity variance
    end
    
    methods (Static)
        % generates the covariance of process noise matrix
        function Q = Qgen(var,dt)
            Qp = [[var.P*eye(3), zeros(3,3)];                     % receiver orbit model covariance
                  [zeros(3,3), var.V*eye(3)]];
            Qc = [[var.Sb*dt + 0.5*var.Sf*dt^3, 0.5*var.Sf*dt^2]; % receiver clock model covariance
                  [0.5*var.Sf*dt^2, var.Sf*dt]];
            Q = [[Qp, zeros(6,2)];                                % combined covariance of process noise matrix
                 [zeros(2,6), Qc]];
        end
        
        % generates the state transition matrix
        function Phi = stateTrans(X,dt,proMode)
            % generate two-body state matrix
            r = norm(X(1:3)');
            F = [[0 0 0 1 0 0];[0 0 0 0 1 0];[0 0 0 0 0 1];
                [((-PARAM.mu/r^3)+(3*PARAM.mu*X(1)^2/r^5)) (3*PARAM.mu*X(1)*X(2)/r^5) (3*PARAM.mu*X(1)*X(3)/r^5) 0 0 0];
                [(3*PARAM.mu*X(1)*X(2)/r^5) ((-PARAM.mu/r^3)+(3*PARAM.mu*X(2)^2/r^5)) (3*PARAM.mu*X(2)*X(3)/r^5) 0 0 0];
                [(3*PARAM.mu*X(1)*X(3)/r^5) (3*PARAM.mu*X(2)*X(3)/r^5) ((-PARAM.mu/r^3)+(3*PARAM.mu*X(3)^2/r^5)) 0 0 0]];
            if proMode == "TBJ2"
                % add J2 effects to state matrix
                gammaJ2 = (3/2)*PARAM.J2*PARAM.mu*PARAM.rearth^2;
                rx = (X(1)^2+X(2)^2+X(3)^2)^(-5/2);
                z1x = 5*X(3)^2;
                z2x = (rx)^(7/5);
                zx = z1x*z2x;
                B0 = (7*z1x*z2x^(9/7))-(5*z2x);
                B1 = -X(1)*B0;
                B2 = -X(2)*B0;
                B3 = X(3)*((10*z2x)-B0);
                B4 = -10*X(1)*z2x;
                B5 = -10*X(2)*z2x;
                B6 = -10*X(3)*z2x;
                df41 = gammaJ2*(zx-rx+(X(1)*B1));
                df42 = gammaJ2*X(1)*B2;
                df43 = gammaJ2*X(1)*B3;
                df51 = gammaJ2*X(2)*B1;
                df52 = gammaJ2*(zx-rx+(X(2)*B2));
                df53 = gammaJ2*X(2)*B3;
                df61 = gammaJ2*X(3)*(B1-B4);
                df62 = gammaJ2*X(3)*(B2-B5);
                df63 = gammaJ2*(zx-rx+(X(3)*(B3-B6)));
                FJ2 = [[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];
                       [df41 df42 df43 0 0 0];
                       [df51 df52 df53 0 0 0];
                       [df61 df62 df63 0 0 0]];
                F = F+FJ2;
            end
            % generate state transition matrix
            Phi = eye(6);
            for i = 1:1:2
                Phi = Phi + (F.^i)*dt/factorial(i);
            end
            PhiC = [[1 dt];[0 1]];
            Phi = [[Phi zeros(6,2)];[zeros(2,6) PhiC]];
        end
        
        % ...
        function X = cow(ode_function,Xprev,Phi,dt,proMode)
            % propagate position and velocity using Cowell's method (Runge-Kutta method #4)
            h = dt; % stopped dividing by 10 here
            b = [0 0 0; 1/2 0 0; 0 1/2 0; 0 0 1];
            c = [1/6 1/3 1/3 1/6];
            t = 0;
            X = Xprev(1:6);
            while t < dt
                Xi = X;
                for ii = 1:4
                    x_inner = Xi;
                    for j = 1:ii-1
                        x_inner = x_inner + h*b(ii,j)*f(:,j);
                    end
                    f(:,ii) = feval(ode_function,x_inner,proMode); % calls method: dfdt = EKF.rates(X,proMode)
                end
                h = min(h, dt-t);
                t = t + h;
                X = Xi + h*f*c';
            end
            % propagate clock states using state transition matrix
            Xclock = Phi(7:8,7:8)*Xprev(7:8);
            X = [X;Xclock];
        end
        
        % ...
        function dfdt = rates(X,proMode)
            R = X(1:3)';         %     [m], position vector
            V = X(4:6)';         %   [m/s], velocity vector
            r = norm(R);         %     [m], distance from Earth's center
            a = -PARAM.mu*R/r^3; % [m/s^2], two-body acceleration
            if proMode == "TBJ2"
                % add J2 acceleration contribution
                gammaJ2 = (3/2)*PARAM.J2*PARAM.mu*PARAM.rearth^2; % J2-Pert Multiplication Constant [...]
                rx = (X(1)^2+X(2)^2+X(3)^2)^(-5/2);               % Algorithm Variable              [...]
                z1x = 5*X(3)^2;                                   % Algorithm Variable              [...]
                z2x = (rx)^(7/5);                                 % Algorithm Variable              [...]
                zx = z1x*z2x;                                     % Algorithm Variable              [...]
                ax_J2 = gammaJ2*X(1)*(zx-rx);                     % J2 Acceleration: x-component    [...]
                ay_J2 = gammaJ2*X(2)*(zx-rx);                     % J2 Acceleration: y-component    [...]
                az_J2 = gammaJ2*X(3)*(zx-3*rx);                   % J2 Acceleration: z-component    [...]
                aJ2 = [ax_J2 ay_J2 az_J2];                        % J2 Acceleration Vector          [m/s^2]
                a = a + aJ2;
            end
            dfdt = [V a]'; % [m/s]/[m/s^2], differential state vector
        end
        
        % ...
        function [Zr,G,R] = ZGRgen(Z,X,mesMode,mesVar)
            if mesMode == "PR"
                norms = zeros(Z.numSV,1);
                prEst = zeros(Z.numSV,1);
                G = zeros(Z.numSV,8);
                R = mesVar.Pr*eye(Z.numSV);
                for ii = 1:1:Z.numSV
                    norms(ii) = norm(Z.xs(ii,:)-X(1:3)');                          % [m], euclidean distance bt. SV and estimated position
                    prEst(ii) = norms(ii) + X(7);                                  % [m], estimated pseudorange, GPS textbook (ch. 9, eq. 2)
                    G(ii,:) = [-(Z.xs(ii,:)-X(1:3)')/norms(ii), zeros(1,3), 1, 0]; % measurement connection matrix, GPS textbook (pg. 421)
                end
                Zr = Z.psr - prEst; % [m], measurement residual, GPS textbook (pg. 422)
                
            elseif mesMode == "SPS"
                Zr = Z.X - X;                         % calculate measurement residual
                G = eye(6,6);                         % ...
                R = [[mesVar.P*eye(3) zeros(3,3)]; ...
                     [zeros(3,3) mesVar.V*eye(3)]];    % ... 
            end
        end
    end
    
    % ...
    methods (Static)
        % ...
        function [X,P] = statePred(ekf,X,P)
            Phi = EKF.stateTrans(X,ekf.dt,ekf.proMode);       % step 0: construct state transition matrix
            X = EKF.cow(@EKF.rates,X,Phi,ekf.dt,ekf.proMode); % step 1: propagate state matrix
            P = Phi*P*Phi' + ekf.Q;                           % step 2: propagate state covariance matrix
        end
        
        % ...
        function [X,P,Zr] = mesUpdate(ekf,X,P,Z)
            [Zr,G,R] = EKF.ZGRgen(Z,X,ekf.mesMode,ekf.mesVar); % step 3: generate observation variables
            K = P*G'*inv(G*P*G'+R);                            % step 4: compute the gain matrix
            X = X + K*(Zr);                                    % step 5: update the state
            P = (eye(8)-K*G)*P*(eye(8)-K*G)' + K*R*K';         % step 6: update the state covariance
        end
    end
end