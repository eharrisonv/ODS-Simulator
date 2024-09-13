%       title: GPS.m
%      author: Elijah Vautour
%        date: November 6, 2022
% description: The GPS class creates a simulated GPS sensor object which in
%              turn generates GPS observables and single point solution (SPS) state
%              estimates.

classdef GPS < matlab.mixin.Copyable
    % GPS sensor properties
    properties
        ID
        genDate
        data = struct('obsFile',[], ...
                      'obsData',[], ...
                      'ephFile',[], ...
                      'ephData',[])
        prop = struct('simTimeDelta',[], ...
                      'compTime',[])
        obs = struct('T',[], ...
                     'time',[], ...
                     'numSV',[], ...
                     'psr',[], ...
                     'xs',[])              
    end
    
    % GPS sensor equations
    methods (Static)
        % ...
        function utcVec = utcVecGen(tUNIX)
            gpsLeap = GPS.leapCalc(tUNIX);
            tUTC = tUNIX - gpsLeap;
            dateTime = datetime(tUTC,'ConvertFrom','posixTime','TimeZone','UTC','Format','dd-MMM-yyyy HH:mm:ss.SSS');
            utcVec = [dateTime.Year ...
                      dateTime.Month ...
                      dateTime.Day ...
                      dateTime.Hour ...
                      dateTime.Minute ...
                      dateTime.Second];
        end
        
        % calculates necessary time parameters for GPS algorithm
        function time = timeCalc(unixTime)
            time = DATA.timeStruct;
            time.tUNIX = unixTime;
            time.gpsLeap = GPS.leapCalc(unixTime);
            time.tUTC = unixTime - time.gpsLeap;
            time.tGPS = unixTime - PARAM.gps2unix;
            time.gpsWeek = floor(time.tGPS/PARAM.secPerWeek);
            gpsDay = floor(time.tGPS/PARAM.secPerDay);
            time.gpsTow = time.tGPS - (time.gpsWeek*PARAM.secPerWeek);
            time.gpsTod = time.tGPS - (gpsDay*PARAM.secPerDay);
            time.dateTime = datetime(time.tUTC,'ConvertFrom','posixTime','TimeZone','UTC','Format','dd-MMM-yyyy HH:mm:ss.SSS');
            time.utcVec = [time.dateTime.Year ...
                           time.dateTime.Month ...
                           time.dateTime.Day ...
                           time.dateTime.Hour ...
                           time.dateTime.Minute ...
                           time.dateTime.Second];
        end
        
        % determines current number of leap seconds past GPS time epoch
        function gpsLeap = leapCalc(unixTime)
            gpsLeap = 0;
            for i = 1:1:length(PARAM.leapEpoch)
                if (unixTime - PARAM.leapEpoch(i)) < 0
                    break;
                else
                    gpsLeap = i;
                end
            end
        end
        
        % calculates the GPS SV time
        function [ts,tk] = refEpoch(psr,tow,toe)
            tau = psr/PARAM.c; % [s], estimated propagation time
            ts = tow - tau;    % [s], estimated satellite time
            tk = ts - toe;     % [s], estimated time from ephemeris reference epoch (GPS-ICD, pg.##)
            if (tk > 0.5*PARAM.secPerWeek)
                tk = tk - PARAM.secPerWeek;
            elseif (tk < -0.5*PARAM.secPerWeek)
                tk = tk + PARAM.secPerWeek;
            end
        end
        
        % calculates the GPS SV eccentric anomaly orbital parameter
        function Ek = eccAnom(M0,n,tk,ecc)
            Mk = M0 + n*tk;
            syms E;
            eqn = E - ecc*sin(E) == Mk;
            solx = vpasolve(eqn, E);
            Ek = double(solx);
        end
        
        % calculates the GPS SV position and corrects the raw pseudorange measurement
        function [psr,xs] = satPos(eph,towMes,psr)
            % 1. calculate missing ephemeris parameters
            A = eph.sqrtA^2;         %     [m], semi-major axis       (GPS-ICD, pg.##)
            n0 = sqrt(PARAM.mu/A^3); % [rad/s], computed mean motion  (GPS-ICD, pg.##)
            n = n0 + eph.deltaN;     % [rad/s], corrected mean motion (GPS-ICD, pg.##)
            
            % 2. calculate pseudorange SV clock bias correction
            [ts,tk] = GPS.refEpoch(psr,towMes,eph.toe);   % time from ephemeris reference epoch
            Ek = GPS.eccAnom(eph.M0,n,tk,eph.ecc);        % eccentric anomaly
            dtr = PARAM.F*(eph.ecc^eph.sqrtA)*sin(Ek);    % ...
            dts = eph.af0 + eph.af1*(ts-eph.toc) + ...
                  eph.af2*(ts-eph.toc)^2 + dtr - eph.Tgd; % ...
            psr = psr + PARAM.c*dts;                      % pseudorange (SV clock bias correction)
            
            % 3. calculate SV ECEF position
            [~,tk] = GPS.refEpoch(psr,towMes,eph.toe);                 % corrected time from ephemeris reference epoch
            Ek = GPS.eccAnom(eph.M0,n,tk,eph.ecc);                     % corrected eccentric anomaly
            vk = 2*atan(sqrt((1+eph.ecc)/(1-eph.ecc))*tan(Ek/2));      % [rad], true anomaly                          (GPS-ICD, pg.##)
            phi = vk + eph.omega;                                      % [rad], argument of latitude                  (GPS-ICD, pg.##)
            du = eph.Cus*sin(2*phi) + eph.Cuc*cos(2*phi);              % [rad], argument of latitude correction       (GPS-ICD, pg.##)
            dr = eph.Crs*sin(2*phi) + eph.Crc*cos(2*phi);              %   [m], radius correction                     (GPS-ICD, pg.##)
            di = eph.Cis*sin(2*phi) + eph.Cic*cos(2*phi);              % [rad], inclination correction                (GPS-ICD, pg.##)
            u = phi + du;                                              % [rad], corrected argument of latitude        (GPS-ICD, pg.##)
            r = A*(1 - eph.ecc*cos(Ek)) + dr;                          %   [m], corrected radius                      (GPS-ICD, pg.##)
            i = eph.i0 + di + (eph.IDOT)*tk;                           % [rad], corrected inclination                 (GPS-ICD, pg.##)
            xprime = r*cos(u);                                         %   [m], orbital x-position                    (GPS-ICD, pg.##)
            yprime = r*sin(u);                                         %   [m], orbital y-position                    (GPS-ICD, pg.##)
            OMEGAk = eph.OMEGA0 + (eph.raanRate - PARAM.wE)*tk - ...
                     PARAM.wE*eph.toe;                                 % [rad], corrected longitude of ascending node (GPS-ICD, pg.##)
            xk = xprime*cos(OMEGAk) - yprime*cos(i)*sin(OMEGAk);       %   [m], ECEF x-coordinate                     (GPS-ICD, pg.##)
            yk = xprime*sin(OMEGAk) + yprime*cos(i)*cos(OMEGAk);       %   [m], ECEF y-coordinate                     (GPS-ICD, pg.##)
            zk = yprime*sin(i);                                        %   [m], ECEF z-coordinate                     (GPS-ICD, pg.##)
            xs = [xk yk zk];                                           % ...
        end
        
        % rotates the SV position from ECEF coordinates at signal transmission time to
        % ECEF coordinates at signal reception time
        function xsCor = xsCor(xsRaw,prMes)
            xsCor = zeros(length(prMes),3);
            for i = 1:1:length(prMes)
                theta = PARAM.wE*(prMes(i,:)/PARAM.c);
                xsCor(i,:) = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1]*xsRaw(i,:)';
            end
        end
    end
    
    % GPS sensor algorithms
    methods (Static)
        % computes the GPS LS solution
        function [x,b] = kinSol(numSV,prMes,xs,x,b)
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
        end
        
        % GPS sensor update for observation dataset generation
        function [compTime,time,psr,xs] = navGen(msg,t)
            % 1. input raw observables
            %svID = msg{2};
            numSV = msg{1}; % # of tracked GPS space vehicles (SV) or satellites
            prMes = msg{3}; % observable 1: raw pseudorange measurement [m], (numSV x 1)
            eph = msg{4};   % observable 2: SV navigation parameters    [-], {numSV x 1}
            % compute GPS observables
            tic
            time = GPS.timeCalc(t);
            psr = zeros(numSV,1);
            xs = zeros(numSV,3);
            for i = 1:1:numSV
                [psr(i),xs(i,:)] = GPS.satPos(eph{i},time.gpsTow,prMes(i));
            end
            xs = GPS.xsCor(xs,psr);
            compTime = toc;
        end
    end
end