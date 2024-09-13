%       title: DATA.m
%      author: Elijah Vautour
%        date: November 6, 2022
% description: The DATA class contains the functions used for communication
%              and data handling in the OPG.m simulator.

classdef DATA
    % navigation data structures
    properties (Hidden, Constant)
        % GPS time structure
        timeStruct = struct('gpsLeap',[], ...  %    [s], leap seconds added since standGPSepoch
                            'gpsWeek',[], ...  % [week], GPS week
                            'gpsTow',[], ...   %    [s], GPS time of week
                            'gpsTod',[], ...   %    [s], GPS time of day
                            'tGPS',[], ...     %    [s], seconds elapsed since GPS time epoch
                            'tUNIX',[], ...    %    [s], seconds elapsed since UTC time epoch
                            'tUTC',[], ...     %    [s], seconds elapsed since UTC time epoch (leap seconds included)
                            'dateTime',[], ... %    [-], UTC date time string
                            'utcVec',[]);      %    [-], UTC date time vector
        
        % SPS navigation solution data structure
        XkStruct = struct('ID',[], ...
                          'class',[], ...
                          'gps',[], ...
                          'compTime',[], ...
                          'navFlag',[], ...
                          'navOut',[])

        % EKF navigation solution data structure
        XeStruct = struct('ID',[], ...
                          'class',[], ...
                          'ekf',[], ...
                          'gps',[], ...
                          'compTime',[], ...
                          'navFlag',[], ...
                          'navOut',[])
    end
    
    % simulation methods
    methods (Static)
        % path handling for simulation setup
        function [datPath,inPath,outPath] = simSet(rootPath,simData,simDate)
            datPath = strcat(rootPath,"\datasets\",simData,"\",simDate);
            inPath = strcat(datPath,"\preprocessed data\");
            addpath(inPath);
            curDate = datestr(datetime,'yyyy-mm-dd');
            outPath = strcat(rootPath,"\results\",curDate,"\");
        end
        
        % outputs matlab data structures based on path definition
        function simOut(datOut,datOutID,resFold)
            % create results folder if it doesn't already exist
            if exist(resFold,'dir') == 0
                mkdir(resFold)
            end
            % save output data, if filename exists already prompt user to input new name
            while 1
                fid = strcat(resFold,datOutID,".mat");
                if exist(fid,'file') ~= 0
                    prompt = strcat("output filename '",datOutID,"' exists, enter new output filename: ");
                    datOutID = input(prompt,'s');
                    clc
                else
                    fprintf('saving results...\n')
                    save(fid,'-struct','datOut')
                    clc
                    break
                end
            end
        end
        
        % initializes simTime progress tracker
        function simTime = simProgInit(ID,T)
            simTime.ID = ID;
            simTime.ticStart = tic;
            simTime.ticProg = tic;
            simTime.simDuration = T(end)-T(1);
            simTime.simStart = T(1);
            simTime.delta = 0;
            simTime.simProgPrev = 0;
            fprintf("simulation starting...\n")
        end

        % updates simTime progress tracker
        function simTime = simProg(simTime,t)
            tLim = 2.5;
            %pLim = 1.0; % gps generation
            pLim = 10.0; % navigation simulation
            if toc(simTime.ticProg) > tLim
                simProg = 100*(t-simTime.simStart)/simTime.simDuration;
                if (simProg-simTime.simProgPrev) > pLim
                    clc
                    simTime.simProgPrev = floor(simProg/pLim)*pLim;
                    fprintf(strcat("computing '",simTime.ID,"' navigation solution...\n"))
                    fprintf("        time elapsed = %.2f seconds\n",simTime.delta)
                    fprintf(" simulation progress = %.2f percent\n",simTime.simProgPrev)
                end
                simTime.ticProg = tic;
            end
            simTime.delta = toc(simTime.ticStart);
        end
    end
    
    % GPS dataset methods
    methods (Static)
        % GPS sensor: observable generation
        function [gpsFIX,time,numSV,psr,xs] = gpsObs(data,t)
            t0 = data.T(1);
            tf = data.T(end);
            dt = data.T(2)-data.T(1);
            index = 1+((t-t0)/dt);
            if (t < t0 || ...          % timestep outside of range (before t0)
                t > tf || ...          % timestep outside of range (after tf)
                floor(index) ~= index) % no data available for timestep
                gpsFIX = false;
                time = [];
                numSV = [];
                psr = [];
                xs = [];
            else
                gpsFIX = true;
                time = data.time{index};
                numSV = data.numSV(index);
                psr = data.psr{index};
                xs = data.xs{index};
            end
        end
        
        % GPS sensor: object generation
        function gps = gpsGen(ID,obsFile,ephFile,simLength)
            % load raw GPS observation and navigation datasets
            obsData = load(obsFile,'-mat');
            ephData = load(ephFile,'-mat');
            % store ephemeris data in cell structures organized according to observables
            eph = cell(obsData.N,1);
            for i = 1:1:obsData.N
                t = obsData.time(i);
                numSV = obsData.numSV(i);
                svID = obsData.svID{i};
                eph{i} = cell(numSV,1);
                for j = 1:1:numSV
                    index = ephData.(svID(j)).count;
                    while ephData.(svID(j)).t(index) > t
                        if index == 1
                            break
                        else
                            index = index - 1;
                        end
                    end
                    eph{i}{j,1} = ephData.(svID(j)).data{index};
                end
            end
            obsData.eph = eph;
            % generate GPS SV ephemerides for all observation data timesteps
            simN = round((3600*simLength)/obsData.dt);
            if simN > obsData.N
                simN = obsData.N;
            end
            psr = cell(simN,1);
            xs = cell(simN,1);
            time = cell(simN,1);
            compTime = zeros(simN,1);
            simTime = DATA.simProgInit(ID,obsData.time(1:simN)');
            for i = 1:1:simN
                % generate observation data for given timestep
                t = obsData.time(i);
                msg{1} = obsData.numSV(i);
                msg{2} = obsData.svID{i};
                msg{3} = obsData.prMes{i};
                msg{4} = obsData.eph{i};
                % solve for gps observables
                [compTime(i),time{i},psr{i},xs{i}] = GPS.navGen(msg,t);
                simTime = DATA.simProg(simTime,t);
            end
            % store data into GPS object
            gps = GPS; % create GPS object
            gps.ID = ID;
            gps.genDate = datestr(datetime,'yyyy-mm-dd');
            gps.data.obsFile = obsFile;
            gps.data.obsData = obsData;
            gps.data.ephFile = ephFile;
            gps.data.ephData = ephData;
            gps.prop.simTimeDelta = simTime.delta;
            gps.prop.compTime = compTime;
            gps.obs.T = obsData.time(1:simN)';
            gps.obs.time = time;
            gps.obs.numSV = obsData.numSV(1:simN);
            gps.obs.psr = psr;
            gps.obs.xs = xs;
        end
    end
    
    % NAV dataset methods
    methods (Static)
        % NAV model: truth solution setup
        function Xa = truNavSet(T)
            Xa = load('navData.mat','-mat');
            Xa.class = "truth";
            % ...
            Xa.truNav = datSelect(Xa.truNav);
            Xa.navOut = datSelect(Xa.data.GNV1A.data);
            % ...
            Xa.data.GNV1A = Xa.data.GNV1A.input;
            Xa.data.GNV1B = Xa.data.GNV1B.input;
            Xa.data.CLK1B = Xa.data.CLK1B.input;
            % ...
            function data = datSelect(data)
                t0 = data(1,1);
                tf = data(end,1);
                if T(1) >= t0 && T(end) <= tf
                   N = length(data);
                   dt = data(2,1)-t0;
                   ind0 = 1 + ((T(1)-t0)/dt);
                   indf = 1 + ((T(end)-t0)/dt);
                   data = data(ind0:indf,:);
                else
                    fprintf("test data not within range...\n");
                end
            end
        end
        
        % NAV model: SGP4 solution setup
        function Xt = sgp4NavSet(ID,sgp4,compTime,navOut)
            Xt.ID = ID;
            Xt.class = "SGP4";
            Xt.tle = sgp4.tle;
            Xt.compTime = compTime;
            Xt.navOut = navOut;
            clc % comment out for debugging
        end
        
        % SGP4 model: ephemeris selection function
        function ecef = sgp4dat(data,t)
            index = 1 + ((t - data.t0)/data.dt);
            if (t < data.t0 || ...     % timestep outside of range (before t0)
                t > data.tf || ...     % timestep outside of range (after tf)
                floor(index) ~= index) % no data available for this timestep
                ecef = [];
            else
                ecef = data.ecef(index,:);
            end
        end
        
        % EKF model: object initialization
        function ekf = ekfSet(ekf_settings,proMode,mesMode)
            ekf = EKF;
            ekf_settings_fields = fieldnames(ekf_settings);
            for i = 1:1:length(ekf_settings_fields)
                fieldSet = ekf_settings_fields{i};
                ekf.(fieldSet) = ekf_settings.(fieldSet);        
            end
            ekf.proMode = proMode;
            ekf.mesMode = mesMode;
        end
    end
end