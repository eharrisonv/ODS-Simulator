%       title: OPD.m
%      author: Elijah Vautour
%        date: November 6, 2022
% description: The OPD class contains the functions used to predict flybys
%              with an Earth based targets viewing cone.

classdef OPD
    % OPD_class constant properties
    properties (Hidden, Constant)
        rEarth = 6.3781E6;              %     [m], Earth radius
        wEarth = 7.2921159E-5*(180/pi); % [deg/s], Earth angular velocity
    end
    
    % OPD_class object properties
    properties
        % target definition
        target = struct('lat',[], ...       % [deg], target geodetic latitude
                        'lon',[], ...       % [deg], target geodetic longitude
                        'viewAngle',[], ... % [deg], target horizon half-angle
                        'viewCone',[]);     % [deg], target horizon geodetic coordinates (used for plotting results)

        % algorithm flags
        flag = struct('appDrn',[], ...  % indicates satellite approach direction (1 = ascending, 2 = descending)
                      'predGen',false, ... % indicates satellite flyby prediction generation
                      'newPred',false);    % indicates newly completed satellite flyby predictions (used externally)
                  
        % satellite tracking
        track = struct('geoNow',[], ...   % satellite current geodetic coordinates, [time latitude longitude]
                       'TLL',[0 0], ...   %   [s], most recent TLL crossing times, [t_asc t_dsc]
                       'per',[], ...      %   [s], estimated orbital period, [P_asc P_dsc]
                       'latPrev',0, ...   % [deg], satellite previous geodetic latitude (used to determine approach direction)
                       'latSignPrev',[]); % previous satellite-target relative latitude sign (used to check for TLL crossing events)

        % target prediction tracking
        pred = struct('coarse', ...
                         [999 999], ... % coarse orbit prediction (number of TLL crossings before next two upcoming flybys, ascending and descending)
                      'fine', ...
                         struct('N',[], ...    % orbits remaining until predicted pass
                                'drn',[], ...  % approach direction of predicted pass
                                'rise',[], ... % satellite geodetic coordinates at predicted target horizon rise, [time latitude longitude]
                                'set',[], ...  % satellite geodetic coordinates at predicted target horizon set, [time latitude longitude]
                                'TCSP',[]));   % satellite geodetic coordinates at predicted target-closest-satellite approach, [time latitude longitude]
    end
    
    % OPD_class static methods
    methods (Static)
        % orbital target generation function
        function target = target_gen(tit,lat,lon,R)
            target.tit = tit;
            target.lat = lat;
            target.lon = lon;
            target.viewAngle = acos(OPD.rEarth/R)*180/pi;
            theta = 0:pi/16:2*pi;
            target.viewCone = zeros(length(theta),2);
            for i = 1:1:length(theta)
                target.viewCone(i,1) = target.lat + target.viewAngle*cos(theta(i));
                target.viewCone(i,2) = target.lon + target.viewAngle*sin(theta(i));
            end
        end
    end
    
    % OPD_class object methods
    methods
        % OPD algorithm driving function
        function opd = OPD_drv(opd,geoNow)
            % 1. update current tracking parameters (input block)
            opd.track.geoNow = geoNow;
            latDeltaNow = opd.track.geoNow(2) - opd.target.lat;

            % 2. check satellite approach direction (decision block)
            if opd.track.geoNow(2) > opd.track.latPrev
                opd.flag.appDrn = 1;
            elseif opd.track.geoNow(2) < opd.track.latPrev
                opd.flag.appDrn = 2;
            end

            % 3. check for TLL crossing (decision block)
            if sign(latDeltaNow) ~= opd.track.latSignPrev
                % 3-1. update TLL crossing time and orbital period estimate
                TLL_prev = opd.track.TLL(opd.flag.appDrn);
                if TLL_prev ~= 0 % condition: if previous TLL crossing timestamp available
                    opd.track.per(opd.flag.appDrn) = opd.track.geoNow(1) - TLL_prev;
                end
                opd.track.TLL(opd.flag.appDrn) = opd.track.geoNow(1);
                % 3-2. coarse orbit prediction
                predMode = 2;
                opd.pred.coarse(1,opd.flag.appDrn) = coarsePred(opd,predMode);
            end

            % 4. check for upcoming flyby in next pass (decision block)
            if opd.pred.coarse(opd.flag.appDrn) <= 2 && ... % condition 1: if target flyby predicted on next ascending or descending pass
               abs(latDeltaNow) < opd.target.viewAngle      % condition 2: if satellite is currently within set target latitude bounds
                % 4-1. fine orbit prediction
                N = 1;
                opd = finePred(opd,1);
            end

            % 5. update previous tracking parameters (output block)
            opd.track.latPrev = opd.track.geoNow(2);
            opd.track.latSignPrev = sign(latDeltaNow);
        end
        
        % coarse orbit prediction function
        function N = coarsePred(opd,predMode)
            % ...
            geoNow = opd.track.geoNow;
            per = opd.track.per(opd.flag.appDrn);
            rot = per*opd.wEarth;
            
            % ...
            if predMode == 1
                if geoNow(3) > opd.target.lon + opd.target.viewAngle
                    dvv = geoNow(3) - (opd.target.lon + opd.target.viewAngle);
                else
                    dvv = 360 + geoNow(3) - (opd.target.lon + opd.target.viewAngle);
                end
                N = floor(abs(dvv/rot))+1;
            else
                if geoNow(3) > opd.target.lon + opd.target.viewAngle
                    dvv = geoNow(3) - (opd.target.lon + opd.target.viewAngle);
                elseif geoNow(3) - rot > opd.target.lon - opd.target.viewAngle
                    dvv = 0;
                else
                    dvv = 360 + geoNow(3) - (opd.target.lon + opd.target.viewAngle);
                end
                N = floor(abs(dvv/rot))+1;
            end
        end
        
        % fine orbit prediction function
        function opd = finePred(opd,N)
            % predict geodetic coordinates N orbits in the future
            geoNext = opd.track.geoNow;
            geoNext(1) = geoNext(1) + N*opd.track.per(opd.flag.appDrn);
            geoNext(3) = geoNext(3) - N*(opd.wEarth*opd.track.per(opd.flag.appDrn));
            while geoNext(3) < -180
                geoNext(3) = geoNext(3) + 360;
            end
            latDeltaNext = geoNext(2) - opd.target.lat;
            lonDeltaNext = geoNext(3) - opd.target.lon;
            rangeNext = norm([latDeltaNext lonDeltaNext]);
            % check if predicted coordinates are within range of target
            if rangeNext < opd.target.viewAngle
                % initialize fine prediction estimate (rising horizon pass)
                if opd.flag.predGen == false
                    opd.flag.predGen = true;
                    opd.pred.fine.N = N;
                    opd.pred.fine.drn = opd.flag.appDrn;
                    opd.pred.fine.rise = geoNext;
                    opd.pred.fine.TCSP = geoNext;
                % update fine prediction estimate (TCSP and setting horizon pass)
                else
                    rangeTCSP = norm([opd.pred.fine.TCSP(2) - opd.target.lat, ...
                                      opd.pred.fine.TCSP(3) - opd.target.lon]);
                    if rangeNext < rangeTCSP
                        opd.pred.fine.TCSP = geoNext;
                    end
                end
                opd.pred.fine.set = geoNext;
            else
                % update flags for newly completed prediction
                if opd.flag.predGen == true
                    opd.flag.predGen = false;
                    opd.flag.newPred = true;
                end
            end
        end
    end
    
end