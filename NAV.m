%       title: NAV.m
%      author: Elijah Vautour
%        date: November 6, 2022
% description: The NAV class contains the actual navigation algorithms which are
%              tested in the OPG.m simulator.

classdef NAV
    % common navigation algorithm methods
    methods (Static)
        % computes the GPS LS solution for tracked state variables
        function X = kinSol(numSV,prMes,xs,X,dt)
            x0 = X(1:3)'; x = x0;
            b0 = X(7)'; b = b0;
            dx = 100*ones(1,3);
            while norm(dx) > 1e-3
                norms = sqrt(sum((xs-x).^2,2));     % euclidean distance bt. SV and estimated receiver position
                prEst = norms + b;                  % estimated psuedorange, GPS textbook (ch. 9, eq. 2)
                dp = prEst - prMes;                 % pseudorange residual,  GPS textbook (ch. 9, eq. 3)
                G = [-(xs-x)./norms ones(numSV,1)]; % geometry matrix,       GPS textbook (ch. 9, pg. 413)
                sol = inv(G'*G)*G'*dp;              % estimate correction,   GPS textbook (ch. 9, eq. 9)
                dx = sol(1:3)';
                db = sol(4);
                x = x - dx;
                b = b - db;
            end
            v = (x-x0)/dt;
            f = (b-b0)/dt;
            X = [x v b f]';
        end
        
        % converts state variable from ECEF coordinates to ECI coordinates
        function X = ECEF2ECI(X,tUNIX)
            utcVec = GPS.utcVecGen(tUNIX);
            xECEF = X(1:3);
            vECEF = X(4:6);
            [xECI,vECI] = ecef2eci(utcVec,xECEF,vECEF);
            X = [xECI;vECI;X(7:8)];
        end
        
        % converts state variable from ECI coordinates to ECEF coordinates
        function X = ECI2ECEF(X,tUNIX)
            utcVec = GPS.utcVecGen(tUNIX);
            xECI = X(1:3);
            vECI = X(4:6);
            [xECEF,vECEF] = eci2ecef(utcVec,xECI,vECI);
            X = [xECEF;vECEF;X(7:8)];
        end
        
        % evaluates velocity state convergence
        function validIC = solEval(X,Xprev,dt)
            check = std([X';Xprev']);
            %if mean(check(4:6)) < 100 && dt < 100
            %if mean(check(4:6)) < 300 && dt < 100
            if mean(check(4:6)) < 5000 && dt < 100
                validIC = true;
            else
                validIC = false;
            end
        end
    end
    
    % navigation algorithms
    methods (Static)
        % SGP4 navigation algorithm (currently uses pre-computed ephemeris data)
        function Xt = SGP4(ID,sgp4,T)
            compTime = [];
            navOut = [];
            simTime = DATA.simProgInit(ID,T);
            for t = T
                ecef = DATA.sgp4dat(sgp4,t);
                if ~isempty(ecef)
                    compTime(end+1,:) = 0;
                    navOut(end+1,:) = [t ecef 0 0];
                end
                simTime = DATA.simProg(simTime,t);
            end
            
            % create output data stucture
            Xt = DATA.sgp4NavSet(ID,sgp4,compTime,navOut);
        end
        
        % SGP4 orbit prediction test algorithm
        function Xo = SGP4_Pred(ID,opd,sgp4,T)
            navOut = [];
            geoOut = [];
            geoOutPred = [];
            tPredOut = [0 0 T(1) + 120];
            simTime = DATA.simProgInit(ID,T);
            for t = T
                % ...
                ecef = DATA.sgp4dat(sgp4,t);
                if ~isempty(ecef)
                    % ...
                    Xnow = ecef';
                    llaNow = ecef2lla(Xnow(1:3)','wgs84');
                    geoNow = [t llaNow(1) llaNow(2)];

                    % compute future state estimate
                    opd = OPD_drv(opd,geoNow);
                    geoNext = geoNow;
                    geoNext(1) = geoNow(1) + floor(opd.track.per(opd.flag.appDrn));
                    geoNext(3) = geoNext(3) - (opd.wEarth*floor(opd.track.per(opd.flag.appDrn)));
                    %llaNext = [geoNext(2) geoNext(3)];
                    %Xnext = lla2ecef(llaNext,'wgs84')

                    % ...
                    navOut(end+1,:) = [t Xnow(1:3)'];
                    geoOut(end+1,:) = geoNow;
                    geoOutPred(end+1,:) = geoNext;
                    
                    if opd.flag.newPred == true && ... % condition 1: new prediction available
                       t > tPredOut(3)                    % condition 2: finished previous target flyby
                        opd.flag.newPred = false;
                        tPredOut(end+1,1) = opd.pred.fine.rise(1);
                        tPredOut(end,2) = opd.pred.fine.TCSP(1);
                        tPredOut(end,3) = opd.pred.fine.set(1);
                        %tPredOut(end,:)
                    end
                end
                
                simTime = DATA.simProg(simTime,t);
            end
            
            % create output data structure
            Xo.ID = ID;
            Xo.class = "SGP4";
            Xo.tle = sgp4.tle;
            Xo.target = opd.target;
            Xo.navOut = navOut;
            Xo.geoOut = geoOut;
            Xo.geoOutPred = geoOutPred;
            Xo.tPredOut = tPredOut(2:end,:);
        end
        
        % SPS navigation algorithm
        function Xk = SPS(ID,gps,T)
            % initialize state tracking variables
            tprev = 0;
            Xprev = zeros(8,1);
            % initialize flags
            validIC = false;
            gpsFIX = false;
            solOUT = false;
            % compute navigation solution for given time-series
            gpsFlag = [];
            navOut = [];
            simTime = DATA.simProgInit(ID,T);
            for t = T
                % poll GPS every timestep
                [gpsFIX,~,numSV,psr,xs] = DATA.gpsObs(gps.obs,t);
                % compute LS solution if gpsFIX flag is set
                if gpsFIX == true
                    X = NAV.kinSol(numSV,psr,xs,Xprev,t-tprev);
                    solOUT = true;
                    validIC = NAV.solEval(X,Xprev,t-tprev);
                end
                % store new data if updated solution available
                if solOUT == true
                    gpsFlag(end+1,:) = 1;
                    navOut(end+1,:) = [t X'];
                    tprev = t; Xprev = X; % update apriori states
                    solOUT = false; % reset flag
                end
                simTime = DATA.simProg(simTime,t);
            end
            
            % create output data structure
            Xk.ID = ID;
            Xk.class = "SPS";
            Xk.gps = gps;
            Xk.gpsFlag = gpsFlag;
            Xk.navOut = navOut;
        end
        
        % EKF navigation algorithm
        function Xe = EKF(ID,sgp4,ekf,gps,ICmode,T)
            % initialize EKF variables
            ekf.Q = EKF.Qgen(ekf.proVar,ekf.dt);
            % initialize state tracking variables
            tprev = 0;
            tMesPrev = 0;
            Xprev = zeros(8,1);
            Pprev = eye(8);
            % initialize flags
            validIC = false;
            mesUP = false;
            gpsFIX = false;
            solOUT = false;
            % compute navigation solution for given time-series
            gpsFlag = [];
            mesResOut = [];
            navOut = [];
            simTime = DATA.simProgInit(ID,T);
            for t = T
                % compute navigation solution at specified rate
                testDT = round((t-tprev)*1000)/1000;
                if testDT >= ekf.dt
                    
                    % START OF REAL-TIME NAVIGATION ALGORITHM
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % EKF navigation: IC generation
                    if validIC == false
                        % ICgen_v1: SGP4 solution
                        if ICmode == 1
                            [gpsFIX,time,numSV,psr,xs] = DATA.gpsObs(gps.obs,t);
                            if gpsFIX == true
                                X = [DATA.sgp4dat(sgp4,t) 0 0]';
                                P = Pprev;
                                mesResOut(end+1,:) = [t 0];
                                validIC = true;
                                solOUT = true;
                                mesUP = false;
                            end
                        % ICgen_v2: SPS solution + velocity convergence evaluation
                        elseif ICmode == 2
                            [gpsFIX,time,numSV,psr,xs] = DATA.gpsObs(gps.obs,t);
                            if gpsFIX == true
                                X = NAV.kinSol(numSV,psr,xs,Xprev,t-tprev);
                                P = Pprev;
                                mesResOut(end+1,:) = [t 0];
                                validIC = NAV.solEval(X,Xprev,t-tprev); % evaluate solution
                                solOUT = true;
                                mesUP = true;
                            end
                        % ICgen_v3: SPS/SGP4 combined solution + velocity convergence evaluation
                        elseif ICmode == 3
                            [gpsFIX,time,numSV,psr,xs] = DATA.gpsObs(gps.obs,t);
                            if gpsFIX == true
                                X = NAV.kinSol(numSV,psr,xs,Xprev,t-tprev);
                                P = Pprev;
                                ecef = DATA.sgp4dat(sgp4,t);
                                X(4:6) = ecef(4:6)';
                                mesResOut(end+1,:) = [t 0];
                                validIC = NAV.solEval(X,Xprev,t-tprev); % evaluate solution
                                solOUT = true;
                                mesUP = true;
                            end
                        end
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % EKF navigation: regular operation
                    else
                        % EKF state prediction
                        XprevECI = NAV.ECEF2ECI(Xprev,tprev);      % covert to ECI for state prediction
                        [X,P] = EKF.statePred(ekf,XprevECI,Pprev); % compute state prediction
                        X = NAV.ECI2ECEF(X,t);                     % convert back to ECEF
                        % EKF measurement update
                        if t - tMesPrev > ekf.dtMes
                            [gpsFIX,time,numSV,psr,xs] = DATA.gpsObs(gps.obs,t);
                            % compute measurement update if gpsFIX flag is set
                            if gpsFIX == true
                                Z.psr = psr;
                                Z.xs = xs;
                                Z.numSV = length(psr);
                                [X,P,Zr] = EKF.mesUpdate(ekf,X,P,Z); % compute measurement update
                                mesResOut(end+1,:) = [t norm(Zr)];
                                mesUP = true;
                                tMesPrev = t;
                            end
                        end
                        solOUT = true;
                    end
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % store new data if updated solution available
                    if solOUT == true
                        gpsFlag(end+1,:) = mesUP;
                        navOut(end+1,:) = [t X'];
                        tprev = t; Xprev = X; Pprev = P; % update apriori states
                        solOUT = false; mesUP = false;   % reset flags
                    end
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % END OF REAL-TIME NAVIGATION ALGORITHM
                    
                    simTime = DATA.simProg(simTime,t);
                end
            end
            
            % create output data structure
            Xe.ID = ID;
            Xe.class = "EKF";
            Xe.ekf = ekf;
            Xe.gps = gps;
            Xe.gpsFlag = gpsFlag;
            Xe.mesResOut = mesResOut;
            Xe.navOut = navOut;
        end
        
        % OGNS (EKF version)
        function Xo = OGNS_EKF(ID,opd,sgp4,ekf,gps,ICmode,T)
            % initialize OGNS variables
            ekf.Q = EKF.Qgen(ekf.proVar,ekf.dt);
            % initialize state tracking variables
            tprev = 0;
            tMesPrev = 0;
            Xprev = zeros(8,1); X = Xprev; Xout = X;
            Pprev = eye(8); P = Pprev;
            tPred = [0 0 T(1) + 120];
            % initialize flags
            navMode = 0;
            validIC = false;
            mesUP = false;
            gpsFIX = false;
            solOUT = false;
            % compute navigation solution for given time series
            gpsFlag = [];
            navFlag = [];
            navOut = [];
            tPredOut = [0 0 T(1) + 120];
            simTime = DATA.simProgInit(ID,T);
            for t = T
                testDT = round((t-tprev)*1000)/1000;
                if testDT >= ekf.dt
                    
                    % START OF REAL-TIME NAVIGATION ALGORITHM
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % set navigation mode
                    if t > tPred(1) && t < tPred(3)
                        navMode = 1;     % set fine navigation mode (SPS or EKF)
                    else
                        navMode = 0;     % set coarse navigation mode (SGP4)
                        validIC = false; % reset IC flag
                    end
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % compute state estimate (coarse orbit determination)
                    if navMode == 0
                        % ...
                        ecef = DATA.sgp4dat(sgp4,t);
                        if ~isempty(ecef)
                            X = [ecef'; Xprev(7) + (t - tprev)*Xprev(8); Xprev(8)];
                            Xout = X;
                            solOUT = true;
                        end
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % compute state estimate (fine orbit determination)
                    elseif navMode == 1
                        % EKF navigation: IC generation
                        if validIC == false
                            % ICgen_v1: SPS solution + velocity convergence evaluation
                            if ICmode == 1
                                [gpsFIX,time,numSV,psr,xs] = DATA.gpsObs(gps.obs,t);
                                if gpsFIX == true
                                    X = NAV.kinSol(numSV,psr,xs,Xprev,t-tprev);
                                    P = Pprev;
                                    validIC = NAV.solEval(X,Xprev,t-tprev); % evaluate solution
                                    solOUT = true;
                                    mesUP = true;
                                end
                            % ICgen_v2: SPS/SGP4 combined solution + velocity convergence evaluation
                            else
                                [gpsFIX,time,numSV,psr,xs] = DATA.gpsObs(gps.obs,t);
                                if gpsFIX == true
                                    X = NAV.kinSol(numSV,psr,xs,Xprev,t-tprev);
                                    P = Pprev;
                                    ecef = DATA.sgp4dat(sgp4,t);
                                    X(4:6) = ecef(4:6)';
                                    validIC = NAV.solEval(X,Xprev,t-tprev); % evaluate solution
                                    Xout = X;
                                    solOUT = true;
                                    mesUP = true;
                                end
                            end
                        % EKF navigation: regular operation
                        else
                            % EKF state prediction
                            XprevECI = NAV.ECEF2ECI(Xprev,tprev);      % covert to ECI for state prediction
                            [X,P] = EKF.statePred(ekf,XprevECI,Pprev); % compute state prediction
                            X = NAV.ECI2ECEF(X,t);                     % convert back to ECEF
                            Zr = 0;
                            % EKF measurement update
                            if t - tMesPrev > ekf.dtMes
                                [gpsFIX,time,numSV,psr,xs] = DATA.gpsObs(gps.obs,t);
                                % compute measurement update if gpsFIX flag is set
                                if gpsFIX == true
                                    Z.psr = psr;
                                    Z.xs = xs;
                                    Z.numSV = length(psr);
                                    [X,P,Zr] = EKF.mesUpdate(ekf,X,P,Z); % compute measurement update
                                    tMesPrev = t;
                                end
                            end
                            % EKF solution evaluation
                            validIC = NAV.solEval(X,Xprev,t-tprev); % evaluate solution
                            Xout = X;
                            solOUT = true;
                        end
                    end
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % OGNS guidance: orbit prediction
                    lla = ecef2lla(X(1:3)','wgs84');
                    geoNow = [t lla(1) lla(2)];
                    opd = OPD_drv(opd,geoNow);
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % ...
                    if opd.flag.newPred == true && ... % condition 1: new prediction available
                       t > tPred(3)                    % condition 2: finished previous target flyby
                        opd.flag.newPred = false;
                        tPred(1) = opd.pred.fine.rise(1)-60;
                        tPred(2) = opd.pred.fine.TCSP(1);
                        tPred(3) = opd.pred.fine.set(1)+60;
                        tPredOut(end+1,:) = tPred;
                    end
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % ...
                    if solOUT == true
                        gpsFlag(end+1,:) = mesUP;
                        navFlag(end+1,:) = navMode;
                        navOut(end+1,:) = [t Xout'];
                        tprev = t; Xprev = X; Pprev = P; % update apriori states
                        solOUT = false; mesUP = false;   % reset flags
                    end
                    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % END OF REAL-TIME NAVIGATION ALGORITHM
                    
                    simTime = DATA.simProg(simTime,t);
                end
            end
            
            % create output data structure
            Xo.ID = ID;
            Xo.class = "OGNS";
            Xo.tle = sgp4.tle;
            Xo.target = opd.target;
            Xo.ekf = ekf;
            Xo.gps = gps;
            Xo.gpsFlag = gpsFlag;
            Xo.navFlag = navFlag;
            Xo.navOut = navOut;
            Xo.tPredOut = tPredOut(2:end,:);
        end
        
        % DYN navigation algorithm
        function Xd = DYN(ID,ekf,T,X0)
            % ...
            tprev = T(1);
            Xprev = X0;
            % ...
            navOut = [tprev Xprev'];
            simTime = DATA.simProgInit(ID,T);
            for t = T
                if (t-tprev) >= ekf.dt
                    % ...
                    XprevECI = NAV.ECEF2ECI(Xprev,tprev);
                    Phi = EKF.stateTrans(XprevECI,ekf.dt,ekf.proMode);
                    X = EKF.cow(@EKF.rates,XprevECI,Phi,ekf.dt,ekf.proMode);
                    X = NAV.ECI2ECEF(X,t);
                    % ...
                    gpsFlag(end+1,:) = 0;
                    navOut(end+1,:) = [t X'];
                    tprev = t; Xprev = X;
                end
                simTime = DATA.simProg(simTime,t);
            end
            % ...
            Xd.ID = ID;
            Xd.class = "DYN";
            Xd.dt = ekf.dt;
            Xd.proMode = ekf.proMode;
            Xd.navOut = navOut;
        end
        
    end
end