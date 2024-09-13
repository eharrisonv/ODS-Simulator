%       title: PARAM.m
%      author: Elijah Vautour
%        date: November 6, 2022
% description: The PARAM class contains the list of constant parameters
%              used in the OPG.m simulator including natural constants and
%              mathematical conversions.

classdef PARAM
    properties (Constant)
        % natural constants
        c = 2.99792458e8     %     [m/s], NAT: speed of light
        F = -4.442807633e-10 % [s/m^0.5], GPS: time constant
        J2 = 0.00108263      %       [-], EARTH: second zonal harmonic term
        mu = 3.986005e14     % [m^3/s^2], EARTH: gravitational parameter
        rearth = 6.371e6     %       [m], EARTH: average radius
        wE = 7.2921151467e-5 %   [rad/s], EARTH: rotation rate
        
        % ...
        wgs84 = wgs84Ellipsoid('meters')
        
        % time-frame constants
        secPerWeek = 604800
        secPerDay = 86400
        grace2gps = 630763200
        gps2unix = 315964800
        
        % leap second chart (https://www.ietf.org/timezones/data/leap-seconds.list)
        %  -> each instance corresponds to the posix time at the epoch
        %     which defines the next leap second addition
        %      -> e.g. posixtime(datetime('1981-07-01 00:00:00')) = 362793600
        leapEpoch = [362793600;394329600;425865600;489024000; ...
                     567993600;631152000;662688000;709948800; ...
                     741484800;773020800;820454400;867715200; ...
                     915148800;1136073600;1230768000;1341100800; ...
                     1435708800;1483228800]
    end
end
        