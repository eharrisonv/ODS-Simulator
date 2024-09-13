%       title: PARSE Class Definition
%          by: Elijah Vautour
%        date: November 6, 2022
% description: The PARSE class is used to parse the raw GRACE-FO datasets
%              into MATLAB data structures for testing.

classdef PARSE
    % GRACE-FO file properties
    properties (Hidden, Constant)
        % ...
        gpsNav_structure = struct('ID',[], ...
                                  'class',[], ...
                                  'data',struct('GNV1A',struct('input',[], ...
                                                               'data',[]), ...
                                                'GNV1B',struct('input',[], ...
                                                               'data',[]), ...
                                                'CLK1B',struct('input',[], ...
                                                               'data',[])), ...
                                  'truNav',[], ...
                                  'navOut',[])
        
        % ...
        gpsObs_structure = struct('input',[], ...
                                  'N',[], ...
                                  't0',[], ...
                                  'tf',[], ...
                                  'dt',[], ...
                                  'time',[], ...
                                  'numSV',[], ...
                                  'svID',[], ...
                                  'prMes',[], ...
                                  'cpMes',[])
    end
    
    % GRACE-FO file parsing methods
    methods (Static)
        % convert GRACE-FO epoch time to UNIX epoch time
        function unixTime = grace2unix(graceTime)
            graceGPSepoch = 946728000;            % GRACE-FO GPS time defined by seconds past January 1, 2000, 12:00:00
            unixTime = graceTime + graceGPSepoch; % calculate unixTime, not including leap seconds
        end
        
        % ...
        function [data,input] = readGrace(datPath,datID)
            % ...
            inPath = strcat(datPath,"\raw data");
            fList = ls(inPath);
            for i = 1:1:length(fList)
                if contains(fList(i,1:length(datID)),datID)
                    input = strtrim(fList(i,:));
                    break 
                end
            end
            fid = fopen(strcat(inPath,"\",input));
            % make sure the file name was valid
            if fid == -1
                errordlg(['The file ''' input ''' does not exist.']);
                return
            end
            % skip through header
            end_of_header = 0;
            while end_of_header == 0
                current_line = fgetl(fid);
                if contains(current_line,'# End of YAML header')
                    end_of_header = 1;
                end
            end
            % ...
            data = {};
            while feof(fid)~= 1
                current_line = string(fgetl(fid));
                data{end+1,1}(1,:) = split(current_line);
            end
            fclose(fid);
        end
        
        % ...
        function graceObs(datPath,datID)
            fprintf("parsing GRACE-FO observation data...\n")
            % read raw GRACE observation data
            [datIn,input] = PARSE.readGrace(datPath,datID);
            % parse relevent data
            data = cell(1,4);
            for i = 1:1:length(datIn)
                dataIn_ = datIn{i};
                % check flag for good data
                if dataIn_(6) == "0000111111111111"
                    data{1}(end+1,:) = PARSE.grace2unix(double(dataIn_(1)));        % tGPS
                    data{2}(end+1,:) = string(num2str(double(dataIn_(4)),'G%02u')); % SV PRN
                    data{3}(end+1,:) = double(dataIn_(8));                          % prMesCA
                    data{4}(end+1,:) = double(dataIn_(11));                         % cpMesCA
                end
            end
            % create obsData variable with timestamped GPS observables
            obsData = PARSE.gpsObs_structure;
            obsData.input = input;
            [obsData.time,rawIndex] = unique(data{1});
            obsData.N = length(obsData.time);
            obsData.t0 = obsData.time(1);
            obsData.tf = obsData.time(end);
            obsData.dt = (obsData.tf-obsData.t0)/(obsData.N-1);
            obsData.numSV = zeros(obsData.N,1);
            obsData.svID = cell(obsData.N,1);
            obsData.prMes = cell(obsData.N,1);
            obsData.cpMes = cell(obsData.N,1);
            for i = 1:1:obsData.N
                if i < obsData.N
                    index = rawIndex(i):rawIndex(i+1)-1;
                else
                    index = rawIndex(i):length(data{2});
                end
                obsData.numSV(i) = length(index);
                obsData.svID{i}(:,1) = data{2}(index);
                obsData.prMes{i}(:,1) = data{3}(index);
                obsData.cpMes{i}(:,1) = data{4}(index);
            end
            % output obsData variable
            outPath = strcat(datPath,"\preprocessed data\obsData_",datID,".mat");
            save(outPath,'-struct','obsData')
            clc
        end
        
        % ...
        function graceNav(datPath)
            fprintf("parsing GRACE-FO navigation data...\n")
            % ...
            Xa = PARSE.gpsNav_structure;
            Xa.ID = "GRACE-FO ground-truth";
            % ...
            [datIn,input] = PARSE.readGrace(datPath,'GNV1A');
            navOut = [];
            for i = 1:1:length(datIn)
                datIn_ = datIn{i};
                navOut(i,1) = PARSE.grace2unix(double(datIn_(1)));
                navOut(i,2:4) = double(datIn_(7:9));
                navOut(i,5:7) = double(datIn_(13:15));
                navOut(i,8) = PARAM.c*double(datIn_(19));
                navOut(i,9) = PARAM.c*double(datIn_(21));
            end
            Xa.data.GNV1A.input = input;
            Xa.data.GNV1A.data = navOut;
            % ...
            [datIn,input] = PARSE.readGrace(datPath,'GNV1B');
            truEcef = [];
            for i = 1:1:length(datIn)
                datIn_ = datIn{i};
                truEcef(i,1) = PARSE.grace2unix(double(datIn_(1)));
                truEcef(i,2:4) = double(datIn_(4:6));
                truEcef(i,5:7) = double(datIn_(10:12));
            end
            Xa.data.GNV1B.input = input;
            Xa.data.GNV1B.data = truEcef;
            % ...
            [datIn,input] = PARSE.readGrace(datPath,'CLK1B');
            truClock = [];
            for i = 1:1:length(datIn)
                datIn_ = datIn{i};
                truClock(i,1) = PARSE.grace2unix(double(datIn_(1)));
                truClock(i,2) = double(datIn_(4));
                truClock(i,3) = double(datIn_(6));
            end
            Xa.data.CLK1B.input = input;
            Xa.data.CLK1B.data = truClock;
            % ...
            Xa.truNav = truNavGen(truEcef,truClock);
            outPath = strcat(datPath,"\preprocessed data\navData.mat");
            save(outPath,'-struct','Xa')
            clc
            function truNav = truNavGen(truEcef,truClock)
                % step 0: make sure truCLOCK data is unique
                [~,index] = unique(truClock(:,1));
                truClock = truClock(index,:);
                % step 1: convert CLK1B data timestamps from RCV time to GPS time
                Xq(:,1) = truClock(:,1) + truClock(:,2); % tGPS = tRCV + eps_time;
                Xq(:,2:3) = -PARAM.c*truClock(:,2:3);    % b = -c*eps_time, f = -c*eps_time_drift
                % step 2: interpolate CLK1B data to match GNV1B timestamps
                %  -> linear interpolation of clock bias and drift
                Yq(:,1) = truEcef(:,1); % desired time array
                Yq(:,2:3) = interpn(Xq(:,1),Xq(:,2:3),Yq(:,1),'linear');
                % step 3: convert GNV1B timestamps from GPS time to RCV time
                Xq2 = truEcef;
                Xq2(:,1) = Yq(:,1) + Yq(:,2)/PARAM.c; % tRCV = tGPS - eps_time = tGPS + b/c
                % step 4: interpolate GNV1B data to original timestamps in new format
                %  -> spline interpolation of ECEF position and velocity
                Yq2(:,1) = truEcef(:,1);
                Yq2(:,2:7) = interpn(Xq2(:,1),Xq2(:,2:7),Yq2(:,1),'spline');
                % step 5: trim original CLK1B data to match GNV1B timesteps
                Xq3 = truClock(:,1);
                Xq3(:,2:3) = -PARAM.c*truClock(:,2:3);
                Yq3(:,1:2) = interpn(Xq3(:,1),Xq3(:,2:3),Yq2(:,1),'linear');
                % step 6: create ground-truth navigation solution in RCV time
                truNav = [Yq2 Yq3];
            end
        end
        
        % ...
        function graceSGP4(datPath,simDate)
            fprintf("generating GRACE-FO SGP4 ephemeris...\n")
            % ...
            dt = 1;  % default
            %dt = 10; % for 2018-08-21 to 2018-08-27 (quicker generation)
            t0_UNIX = posixtime(datetime(simDate,'InputFormat','yyyy-MM-dd')); % initial UNIX time timestamp
            t0 = t0_UNIX - GPS.leapCalc(t0_UNIX); % initial UTC time timestamp
            %tf = t0 + 600; % small ephemeris generation for testing
            tf = t0 + 86400 - 1;
            startTime = datetime(t0,'ConvertFrom','posixTime','TimeZone','UTC','Format','dd-MMM-yyyy HH:mm:ss.SSS');
            stopTime = datetime(tf,'ConvertFrom','posixTime','TimeZone','UTC','Format','dd-MMM-yyyy HH:mm:ss.SSS');
            % ...
            datID = 'graceTLE';
            inPath = strcat(datPath,"\raw data");
            fList = ls(inPath);
            for i = 1:1:height(fList)
                if contains(fList(i,1:length(datID)),datID)
                    tleFile = string(strtrim(fList(i,:)));
                    tlePath = strcat(inPath,"\",tleFile);
                    break 
                end
            end
            % generate satellite ephemeris: ECI coordinates in UTC time
            sc = satelliteScenario(startTime,stopTime,dt);                               % create satellite scenario
            graceSAT = satellite(sc,tlePath,"Name","graceSAT","OrbitPropagator","sgp4"); % create a satellite object
            [xECI,vECI,dateTime] = states(graceSAT); % get satellite ECI state vectors
            tUTC = posixtime(dateTime)';             % get UTC timestamps for satellite ephemeris
            % generate satellite ephemeris: ECEF coordinates in UNIX time
            N = length(dateTime);
            xECEF = zeros(size(xECI));
            vECEF = zeros(size(vECI));
            for i = 1:1:N
                utcVec = utcVecGen(dateTime(i));
                [xECEF(:,i),vECEF(:,i)] = eci2ecef(utcVec,xECI(:,i),vECI(:,i));
            end
            T = tUTC + GPS.leapCalc(t0_UNIX);
            Xt.tle = tleData(tleFile,tlePath);
            Xt.t0 = T(1);
            Xt.tf = T(end);
            Xt.dt = dt;
            Xt.T = T;
            Xt.ecef = [xECEF' vECEF'];
            outPath = strcat(datPath,"\preprocessed data\sgp4eph.mat");
            save(outPath,'-struct','Xt')
            clc
            % ...
            function tle = tleData(tleFile,tlePath)
                fid = fopen(tlePath);
                tle.file = tleFile;
                % ...
                line0 = fgetl(fid);
                tle.satID = string(line0); % satellite name identifier
                % ...
                line1 = fgetl(fid);
                year = str2double(line1(19:20));  % UTC element set epoch: last two digits of epoch year
                if year < 57
                    year = string(year + 2000);
                else
                    year = string(year + 1900);
                end
                dayInt = str2double(line1(21:23));  % UTC element set epoch: integer portion of epoch day
                dayFrac = str2double(line1(24:32)); % UTC element set epoch: fractional portion of epoch day
                tle.dateTime = datetime(strcat("1-Jan-",year)) + dayInt + dayFrac - 1;
                % ...
                line2 = fgetl(fid);
                tle.inc = str2double(line2(9:16));               %     [deg], mean orbital element: inclination
                tle.raan = str2double(line2(18:25));             %     [deg], mean orbital element: right ascension of ascending node
                tle.ecc = str2double(strcat("0.",line2(27:33))); %       [-], mean orbital element: eccentricity
                tle.aop = str2double(line2(35:42));              %     [deg], mean orbital element: argument of perigee
                tle.mAnom = str2double(line2(44:51));            %     [deg], mean orbital element: mean anomaly
                tle.mMot = str2double(line2(53:63));             % [rev/day], mean orbital element: mean motion
            end
            % ...
            function utcVec = utcVecGen(dateTime)
                utcVec = [dateTime.Year ...
                          dateTime.Month ...
                          dateTime.Day ...
                          dateTime.Hour ...
                          dateTime.Minute ...
                          dateTime.Second];
            end
        end
    end
    
    % IGS daily GPS broadcast ephemeris file properties
    properties (Hidden, Constant)
        % ephemeris data structure
        ephStruct = struct('weekNo',[], ...   %  [week], GPS week number
                           'Tgd',[], ...      %     [s], group delay differential
                           'af0',[], ...      %     [s], SV clock bias
                           'af1',[], ...      %   [s/s], SV clock drift
                           'af2',[], ...      % [s/s^2], SV clock drift rat
                           'toc',[], ...      %     [s], time of clock
                           'sqrtA',[], ...    % [m^0.5], square root of the semi-major axis
                           'deltaN',[], ...   % [rad/s], mean motion difference from computed value
                           'ecc',[], ...      %     [-], eccentricity
                           'M0',[], ...       %   [rad], mean anomaly at reference time
                           'toe',[], ...      %     [s], time of ephemeris
                           'Crs',[], ...      %     [m], amplitude of the sine correction term to the orbit radius
                           'Cuc',[], ...      %   [rad], amplitude of cosine harmonic correction term to the argument of latitude
                           'Cus',[], ...      %   [rad], amplitude of sine harmonic correction term to the argument of latitude
                           'omega',[], ...    %   [rad], argument of perigee
                           'raanRate',[], ... % [rad/s], rate of right ascension
                           'OMEGA0',[], ...   %   [rad], longitude of ascending node of orbit plane at weekly epoch
                           'i0',[], ...       %   [rad], inclination angle at reference time
                           'IDOT',[], ...     % [rad/s], rate of inclination angle
                           'Cic',[], ...      %   [rad], amplitude of the cosine harmonic correction term to the angle of inclination
                           'Cis',[], ...      %   [rad], amplitude of the sine harmonic correction term to the angle of inclination
                           'Crc',[]);         %     [m], amplitude of the cosine harmonic correction term to the orbit radius
    end
    
    % IGS daily GPS broadcast ephemeris file parsing methods
    methods (Static)
        % open input file, read header
        function [fid,eph] = EPH_fIN(inPath,datID)
            % ...
            fList = ls(inPath);
            for i = 1:1:length(fList)
                if contains(fList(i,1:length(datID)),datID)
                    eph.input = strtrim(fList(i,:));
                    break 
                end
            end
            fid = fopen(strcat(inPath,"\",eph.input));
            % make sure the file name was valid
            if fid == -1
                errordlg(['The file ''' eph.input ''' does not exist.']);
                return;
            end
            % read through header, store IONO parameters
            end_of_header = 0;
            while end_of_header == 0
                current_line = fgetl(fid);
                if contains(current_line,'GPSA')
                    eph.IONO.a0 = double(strcat(current_line(7:13),"E",current_line(15:17)));
                    eph.IONO.a1 = double(strcat(current_line(19:25),"E",current_line(27:29)));
                    eph.IONO.a2 = double(strcat(current_line(31:37),"E",current_line(39:41)));
                    eph.IONO.a3 = double(strcat(current_line(43:49),"E",current_line(51:53)));
                elseif contains(current_line,'GPSB')
                    eph.IONO.b0 = double(strcat(current_line(7:13),"E",current_line(15:17)));
                    eph.IONO.b1 = double(strcat(current_line(19:25),"E",current_line(27:29)));
                    eph.IONO.b2 = double(strcat(current_line(31:37),"E",current_line(39:41)));
                    eph.IONO.b3 = double(strcat(current_line(43:49),"E",current_line(51:53)));
                elseif contains(current_line,'END OF HEADER')
                    end_of_header = 1;
                end
            end
        end
        
        % parse RINEX navigation file
        function ephemeris(datPath)
            fprintf("parsing RINEX navigation data...\n")
            % ...
            inPath = strcat(datPath,"\raw data");
            [fid,ephData] = PARSE.EPH_fIN(inPath,'BRDC00IGS');
            % ...
            while feof(fid) ~= 1
                % parse navigation information
                current_line = fgetl(fid);
                if current_line(1) == 'G' % only use GPS navigation data
                    % identify SV PRN
                    PRN = string(current_line(1:3));
                    % read epoch time
                    Y = double(string(current_line(5:8)));     %  [year], year
                    M = double(string(current_line(10:11)));   % [month], month
                    D = double(string(current_line(13:14)));   %   [day], day
                    H = double(string(current_line(16:17)));   %  [hour], hour
                    min = double(string(current_line(19:20))); %   [min], minute
                    sec = double(string(current_line(22:23))); %     [s], second
                    % read broadcast ephemeris
                    eph_ = PARSE.ephStruct;
                    eph_.af0 = double(string(current_line(24:42))); %     [s], SV clock bias
                    eph_.af1 = double(string(current_line(43:61))); %   [s/s], SV clock drift
                    eph_.af2 = double(string(current_line(62:80))); % [s/s^2], SV clock drift rate
                    [eph_.IODE,eph_.Crs,eph_.deltaN,eph_.M0] = broadcast(fid,1);        % broadcast orbit line 1
                    [eph_.Cuc,eph_.ecc,eph_.Cus,eph_.sqrtA] = broadcast(fid,2);         % broadcast orbit line 2
                    [eph_.toe,eph_.Cic,eph_.OMEGA0,eph_.Cis] = broadcast(fid,3);        % broadcast orbit line 3
                    [eph_.i0,eph_.Crc,eph_.omega,eph_.raanRate] = broadcast(fid,4);     % broadcast orbit line 4
                    [eph_.IDOT,eph_.L2code,eph_.weekNo,eph_.L2flag] = broadcast(fid,5); % broadcast orbit line 5
                    [eph_.svAcc,eph_.svHealth,eph_.Tgd,eph_.IODC] = broadcast(fid,6);   % broadcast orbit line 6
                    [eph_.toc,eph_.fit,~,~] = broadcast(fid,7);                   % broadcast orbit line 7
                    % format and store data
                    if ismember(PRN,fieldnames(ephData))
                        j = ephData.(PRN).count + 1;
                        ephData.(PRN).count = j;
                    else
                        j = 1;
                        ephData.(PRN).count = j;
                        ephData.(PRN).t = [];
                        ephData.(PRN).data = {};
                    end
                    ephData.(PRN).t(1,j) = posixtime(datetime([Y M D H min sec]));
                    ephData.(PRN).data{1,j} = eph_;
                end
            end
            fclose(fid);
            % ...
            outPath = strcat(datPath,"\preprocessed data\ephData.mat");
            save(outPath,'-struct','ephData')
            clc
            % ...
            function [var1, var2, var3, var4] = broadcast(fid,bln)
                current_line = fgetl(fid);
                var1 = double(string(current_line(5:23)));
                var2 = double(string(current_line(24:42)));
                if bln <= 6
                    var3 = double(string(current_line(43:61)));
                    var4 = double(string(current_line(62:80)));
                else
                    var3 = NaN;
                    var4 = NaN;
                end
            end
        end
    end 
end