%       title: OPG.m
%      author: Elijah Vautour
%        date: November 6, 2022
% description: The OPG.m program is the main drive file used to test the
%              different navigation algorithms contained in the NAV class.
% required toolboxes:
%  -> Mapping Toolbox
%  -> Aerospace Toolbox
%  -> Signal Processing Toolbox

%% dataset selection
clear; clc; close all;
format longG

% simulation settings
settings.simData = "GRACE-FO";
% settings.simDate = "2018-08-14";
% settings.simDate = "2018-08-15";
% settings.simDate = "2018-08-16";
% settings.simDate = "2018-08-17";
% settings.simDate = "2018-08-18";
% settings.simDate = "2018-08-19";
% settings.simDate = "2018-08-20";
% settings.simDate = "2018-08-21";
% settings.simDate = "2018-08-22";
% settings.simDate = "2018-08-23";
% settings.simDate = "2018-08-24";
% settings.simDate = "2018-08-25";
% settings.simDate = "2018-08-26";
% settings.simDate = "2018-08-27";
settings.simDate = "2018-08-14_20";
% settings.simDate = "2018-08-14_27";
[settings.datPath,settings.inPath,settings.resFold] = DATA.simSet(pwd,settings.simData,settings.simDate);

%% dataset pre-processing (do once per new dataset)
clc; close all;

% PARSE.graceObs(settings.datPath,'GPS1A') % GRACE-FO data 1: level-1A GPS observables
% PARSE.graceObs(settings.datPath,'GPS1B') % GRACE-FO data 2: level-1B GPS observables
% PARSE.graceNav(settings.datPath);
% PARSE.ephemeris(settings.datPath);
% PARSE.graceSGP4(settings.datPath,settings.simDate);

%% SGP4 ephemeris selection
clc; close all;

% load default sgp4 ephemeris (generated by PARSE)
sgp4 = load("sgp4eph.mat",'-mat');

% load custom sgp4 ephemeris (generated by sgp4gen for 2018-08-14_27 testing)
% sgp4 = load("sgp4eph_14day14tle.mat",'-mat');
% sgp4 = load("sgp4eph_14day1tle.mat",'-mat');

%% GPS sensor simulation (do once per new GPS object)
clc; close all;

% GPS-1 object
% ID = "GRACE-FO C: GPS1A Data, 2018-08-14";
% ephFile = "ephData.mat";
% obsFile = "obsData_GPS1A.mat";
% gpsFile = "gps1.mat";
% simLength = 24; % [h], duration of available GPS observables
% gps1 = DATA.gpsGen(ID,obsFile,ephFile,simLength);
% save(strcat(settings.inPath,"\",gpsFile),'gps1')
load("gps1.mat",'-mat');

%% EKF model setup (currently used models shown here)
clc; close all;

% EKF model settings
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% default EKF settings
ekfSetV1.dt = 1;           %       [s], propagation timestep
ekfSetV1.dtMes = 10;       %       [s], measurement timestep
ekfSetV1.proVar.P = 50;    %     [m^2], eci-position process variance
ekfSetV1.proVar.V = 5;     % [m^2/s^2], eci-velocity process variance
ekfSetV1.proVar.Sf = 30;   %       [?], white noise spectral amplitude
ekfSetV1.proVar.Sb = 0.01; %       [?], white noise spectral amplitude
ekfSetV1.mesVar.P = 150;   %     [m^2], eci-position measurement variance
ekfSetV1.mesVar.V = 75;    % [m^2/s^2], eci-velocity measurement variance
ekfSetV1.mesVar.Pr = 150;  %     [m^2], pseudorange  measurement variance
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% EKF model settings
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% default EKF settings
ekfSetV2.dt = 1;                  %       [s], propagation timestep
ekfSetV2.dtMes = 10;              %       [s], measurement timestep
ekfSetV2.proVar.P = 50;           %     [m^2], eci-position process variance
ekfSetV2.proVar.V = 5;            % [m^2/s^2], eci-velocity process variance
ekfSetV2.proVar.Sf = 4.00000E-19; %       [?], white noise spectral amplitude (GPS textbook, pg. 418)
ekfSetV2.proVar.Sb = 1.15791E-18; %       [?], white noise spectral amplitude (GPS textbook, pg. 418)
ekfSetV2.mesVar.P = 150;          %     [m^2], eci-position measurement variance
ekfSetV2.mesVar.V = 75;           % [m^2/s^2], eci-velocity measurement variance
ekfSetV2.mesVar.Pr = 150;         %     [m^2], pseudorange  measurement variance
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% EKF model: object definitions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ekf1 = DATA.ekfSet(ekfSetV1,"TB","PR");
ekf2 = DATA.ekfSet(ekfSetV1,"TBJ2","PR");
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ekf3 = DATA.ekfSet(ekfSetV2,"TB","PR");
ekf4 = DATA.ekfSet(ekfSetV2,"TBJ2","PR");
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% OPD model setup
clc; close all;

% satellite parameters used to generate target viewing cone
sat.per = (1/sgp4.tle.mMot)*24*3600;         % [s], mean orbital period
sat.R = (PARAM.mu*(sat.per/(2*pi))^2)^(1/3); % [m], mean orbital radius

% Halifax, Canada (44.6476N,63.5728S)
opd1 = OPD;
opd1.target = OPD.target_gen("HFX-CA",44.6476,-63.5728,sat.R);
opd1.track.per = [sat.per sat.per];

% Rio de Janeiro, Brazil (22.9068S,43.1729W)
opd2 = OPD;
opd2.target = OPD.target_gen("RIO-BR",-22.9068,-43.1729,sat.R);
opd2.track.per = [sat.per sat.per];

% Mumbai, India (19.0760N,72.8777E)
opd3 = OPD;
opd3.target = OPD.target_gen("MUM-IN",19.0760,72.8777,sat.R);
opd3.track.per = [sat.per sat.per];

% Quito, Ecuador (0.1807S,78.4678W)
opd4 = OPD;
opd4.target = OPD.target_gen("QUI-EC",-0.1807,-78.4678,sat.R);
opd4.track.per = [sat.per sat.per];

% Nuuk, Greenland (64.1743N,51.7373W)
opd5 = OPD;
opd5.target = OPD.target_gen("NUK-GL",64.1743,-51.7373,sat.R);
opd5.track.per = [sat.per sat.per];

% Concordia Station, Antartica (75.0998S,123.3322E)
opd6 = OPD;
opd6.target = OPD.target_gen("CON-AQ",-75.0998,123.3322,sat.R);
opd6.track.per = [sat.per sat.per];

%% navigation algorithm simulation
clc; close all;

% % TEST 1: SGP4 testing
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% testOut = true;
% settings.navID = "3h_SGP4_test_2018-08-14"; dT = 3; dt = 0.01;
% % settings.navID = "24h_SGP4_test_2018-08-14"; dT = 23.99; dt = 0.01;
% % settings.navID = "96h_SGP4_test_2018-08-14_20"; dT = 96; dt = 0.01;
% % settings.navID = "336h_SGP4_test_14day14tle"; dT = 335.99; dt = 10;
% % settings.navID = "336h_SGP4_test_14day1tle"; dT = 335.99; dt = 10;
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% t0 = gps1.obs.T(1);      % [s], simulation start time
% T = t0 + (0:dt:dT*3600); % [s], simulation time series
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Xnav{1} = DATA.truNavSet(T);
% Xnav{end+1} = NAV.SGP4("SGP4: default ephemeris",sgp4,T);
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % TEST 2: SPS testing
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% testOut = true;
% % settings.navID = "24h_SPS_test_2018-08-14"; dT = 24; dt = 0.01;
% settings.navID = "96h_SPS_test_2018-08-14_20"; dT = 96; dt = 0.01;
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% t0 = gps1.obs.T(1);      % [s], simulation start time
% T = t0 + (0:dt:dT*3600); % [s], simulation time series
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Xnav{1} = DATA.truNavSet(T);
% Xnav{end+1} = NAV.SPS("SPS: GPS1A data",gps1,T);
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % TEST 3: RDEKF IC testing
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% testOut = true;
% settings.navID = "0.5h_EKF_IC_test";
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% dT = 0.5;                 % [h], simulation duration
% dt = 0.01;               % [s], simulation timestep
% t0 = gps1.obs.T(1);      % [s], simulation start time
% T = t0 + (0:dt:dT*3600); % [s], simulation time series
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Xnav{1} = DATA.truNavSet(T);
% Xnav{end+1} = NAV.SGP4("SGP4: default ephemeris",sgp4,T);
% Xnav{end+1} = NAV.SPS("SPS: GPS1A data",gps1,T);
% Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes10, ICmode1",sgp4,ekf1,gps1,1,T);
% Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes10, ICmode1",sgp4,ekf2,gps1,1,T);
% Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes10, ICmode2",sgp4,ekf1,gps1,2,T);
% Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes10, ICmode2",sgp4,ekf2,gps1,2,T);
% Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes10, ICmode3",sgp4,ekf1,gps1,3,T);
% Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes10, ICmode3",sgp4,ekf2,gps1,3,T);
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % TEST 4: RDEKF TB model sensitivity analysis testing
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% testOut = true;
% settings.navID = "3h_EKF_TB_model_test";
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% dT = 3;                  % [h], simulation duration
% dt = 0.01;               % [s], simulation timestep
% t0 = gps1.obs.T(1);      % [s], simulation start time
% T = t0 + (0:dt:dT*3600); % [s], simulation time series
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Xnav{1} = DATA.truNavSet(T);
% ekf1.dtMes = 10; Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes10, ICmode3",sgp4,ekf1,gps1,3,T);
% ekf1.dtMes = 20; Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes20, ICmode3",sgp4,ekf1,gps1,3,T);
% ekf1.dtMes = 30; Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes30, ICmode3",sgp4,ekf1,gps1,3,T);
% ekf1.dtMes = 60; Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes60, ICmode3",sgp4,ekf1,gps1,3,T);
% ekf1.dtMes = 90; Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes90, ICmode3",sgp4,ekf1,gps1,3,T);
% ekf1.dtMes = 180; Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes180, ICmode3",sgp4,ekf1,gps1,3,T);
% ekf1.dtMes = 270; Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes270, ICmode3",sgp4,ekf1,gps1,3,T);
% ekf1.dtMes = 360; Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes360, ICmode3",sgp4,ekf1,gps1,3,T);
% ekf1.dtMes = 450; Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes450, ICmode3",sgp4,ekf1,gps1,3,T);
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % TEST 5: RDEKF TBJ2 model sensitivity analysis testing
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% testOut = true;
% settings.navID = "3h_EKF_TBJ2_model_test";
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% dT = 3;                  % [h], simulation duration
% dt = 0.01;               % [s], simulation timestep
% t0 = gps1.obs.T(1);      % [s], simulation start time
% T = t0 + (0:dt:dT*3600); % [s], simulation time series
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Xnav{1} = DATA.truNavSet(T);
% ekf2.dtMes = 10; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes10, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 20; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes20, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 30; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes30, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 60; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes60, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 90; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes90, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 180; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes180, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 270; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes270, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 360; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes360, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 450; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes450, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 540; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes540, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 630; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes630, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 720; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes720, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 810; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes810, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 900; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes900, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 990; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes990, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 1080; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes1080, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 1170; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes1170, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 1260; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes1260, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 1350; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes1350, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 1440; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes1440, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 1530; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes1530, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 1620; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes1620, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 1710; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes1710, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 1800; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes1800, ICmode3",sgp4,ekf2,gps1,3,T);
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % TEST 6: OGNS testing
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% testOut = true;
% % settings.navID = "23h_OGNS_test"; dT = 23;
% % settings.navID = "23h_OGNS_test_2018-08-15"; dT = 23;
% % settings.navID = "96h_OGNS_test_2018-08-14_20"; dT = 96;
% % settings.navID = "167h_OGNS_test_2018-08-14_20"; dT = 167;
% % settings.navID = "24h_OGNS_MUM_test_2018-08-14_20"; dT = 24;
% % settings.navID = "96h_OGNS_HFX_test2_2018-08-14_20"; dT = 96;
% settings.navID = "96h_OGNS_RIO_test_2018-08-14_20"; dT = 96;
% % settings.navID = "96h_OGNS_MUM_test_2018-08-14_20"; dT = 96;
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% dt = 0.01;               % [s], simulation timestep
% t0 = gps1.obs.T(1);      % [s], simulation start time
% T = t0 + (0:dt:dT*3600); % [s], simulation time series
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Xnav{1} = DATA.truNavSet(T);
% % Xnav{end+1} = NAV.OGNS_EKF("OGNS_EKF: TB/PR, dtMes10, HFX_target",opd1,sgp4,ekf1,gps1,2,T);
% Xnav{end+1} = NAV.OGNS_EKF("OGNS_EKF: TB/PR, dtMes10, RIO_target",opd2,sgp4,ekf1,gps1,2,T);
% % Xnav{end+1} = NAV.OGNS_EKF("OGNS_EKF: TB/PR, dtMes10, MUM_target",opd3,sgp4,ekf1,gps1,2,T);
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % TEST 7: final RDEKF testing
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% testOut = true;
% settings.navID = "24h_EKF_test_2018-08-14"; dT = 24; dt = 0.01;
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% t0 = gps1.obs.T(1);      % [s], simulation start time
% T = t0 + (0:dt:dT*3600); % [s], simulation time series
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Xnav{1} = DATA.truNavSet(T);
% ekf1.dtMes = 10; Xnav{end+1} = NAV.EKF("EKF: TB/PR GPS1A data, dtMes10, ICmode3",sgp4,ekf1,gps1,3,T);
% ekf2.dtMes = 10; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes10, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 810; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes810, ICmode3",sgp4,ekf2,gps1,3,T);
% ekf2.dtMes = 1800; Xnav{end+1} = NAV.EKF("EKF: TBJ2/PR GPS1A data, dtMes1800, ICmode3",sgp4,ekf2,gps1,3,T);
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % TEST 2.1: SGP4 testing
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% testOut = true;
% % settings.navID = "96h_HFX_SGP4_test_2018-08-14_20"; dT = 96; dt = 0.01;
% % settings.navID = "96h_RIO_SGP4_test_2018-08-14_20"; dT = 96; dt = 0.01;
% % settings.navID = "96h_MUM_SGP4_test_2018-08-14_20"; dT = 96; dt = 0.01;
% % settings.navID = "96h_QUI_SGP4_test_2018-08-14_20"; dT = 96; dt = 0.01;
% % settings.navID = "96h_NUK_SGP4_test_2018-08-14_20"; dT = 96; dt = 0.01;
% settings.navID = "96h_CON_SGP4_test_2018-08-14_20"; dT = 96; dt = 0.01;
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% t0 = gps1.obs.T(1);      % [s], simulation start time
% T = t0 + (0:dt:dT*3600); % [s], simulation time series
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Xnav{1} = DATA.truNavSet(T);
% % Xnav{end+1} = NAV.SGP4_Pred("SGP4: default ephemeris",opd1,sgp4,T);
% % Xnav{end+1} = NAV.SGP4_Pred("SGP4: default ephemeris",opd2,sgp4,T);
% % Xnav{end+1} = NAV.SGP4_Pred("SGP4: default ephemeris",opd3,sgp4,T);
% % Xnav{end+1} = NAV.SGP4_Pred("SGP4: default ephemeris",opd4,sgp4,T);
% % Xnav{end+1} = NAV.SGP4_Pred("SGP4: default ephemeris",opd5,sgp4,T);
% Xnav{end+1} = NAV.SGP4_Pred("SGP4: default ephemeris",opd6,sgp4,T);
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% TEST 2.2: OGNS testing
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
testOut = true;
% settings.navID = "96h_OGNS_HFX_test_2018-08-14_20"; dT = 96;
% settings.navID = "96h_OGNS_RIO_test_2018-08-14_20"; dT = 96;
% settings.navID = "96h_OGNS_MUM_test_2018-08-14_20"; dT = 96;
settings.navID = "96h_OGNS_QUI_test_2018-08-14_20"; dT = 96;
% settings.navID = "96h_OGNS_NUK_test_2018-08-14_20"; dT = 96;
% settings.navID = "96h_OGNS_CON_test_2018-08-14_20"; dT = 96;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt = 0.01;               % [s], simulation timestep
t0 = gps1.obs.T(1);      % [s], simulation start time
T = t0 + (0:dt:dT*3600); % [s], simulation time series
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Xnav{1} = DATA.truNavSet(T);
% Xnav{end+1} = NAV.OGNS_EKF("OGNS_EKF: TB/PR, dtMes10, HFX_target",opd1,sgp4,ekf1,gps1,2,T);
% Xnav{end+1} = NAV.OGNS_EKF("OGNS_EKF: TB/PR, dtMes10, RIO_target",opd2,sgp4,ekf1,gps1,2,T);
% Xnav{end+1} = NAV.OGNS_EKF("OGNS_EKF: TB/PR, dtMes10, MUM_target",opd3,sgp4,ekf1,gps1,2,T);
Xnav{end+1} = NAV.OGNS_EKF("OGNS_EKF: TB/PR, dtMes10, QUI_target",opd4,sgp4,ekf1,gps1,2,T);
% Xnav{end+1} = NAV.OGNS_EKF("OGNS_EKF: TB/PR, dtMes10, NUK_target",opd5,sgp4,ekf1,gps1,2,T);
% Xnav{end+1} = NAV.OGNS_EKF("OGNS_EKF: TB/PR, dtMes10, CON_target",opd6,sgp4,ekf1,gps1,2,T);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% output results
if exist('Xnav','var')
    n = length(Xnav);
    IDlist = strings(n,1);
    for i = 1:1:n
        IDlist(i) = Xnav{i}.ID;
    end
    solOut.ID = settings.navID;
    solOut.settings = settings;
    solOut.navList = IDlist;
    solOut.Xnav = Xnav;

    % output solution file
    if testOut == true
        DATA.simOut(solOut,settings.navID,settings.resFold)
    end
end
