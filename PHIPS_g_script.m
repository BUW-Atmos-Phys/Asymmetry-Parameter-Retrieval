function [DATA_out] = PHIPS_g_script(PhipsData, Ng, sizelimit, campaign)
% Calculates g for a group of particles
% Input:    PhipsData  ->  PHIPS level 5 table, data not corrected
%           Ng         ->  Number of particles per group
%           sizelimit  ->  Lower size limit, e.g. 15 micron
%           campaign   ->  Name of the campaign
% Note:     set Ng = 1 if using one averaged ASF
addpath('/Users/emma/Documents/m-codes/PHIPS analysis/AsymmetryFactorAnalysis/src/chebfun-master');


%% Channel correction factor std
if strcmp(campaign,'ACLOUD') || strcmp(campaign,'SOCRATES')
    % Based on individual droplets measured during ARISTO2017
    % from 18 to 42°
    Cm = [1.0313    0.9715    1.3200    1.0924    2.0480    1.3183    1.8443];
    Cs = [0.2656    0.1490    0.2057    0.1718    0.4587    0.3283    0.4339]; 
elseif strcmp(campaign,'CIRRUS-HL')
    % Glass bead calibration results performed during the campaign
    % from 18 to 66°
    Cm = [1.0728    1.0172    1.1056    0.9978    1.0903    1.1346    0.9436    1    1.1575];
    Cs = [0.1185    0.0942    0.1032    0.0965    0.0893    0.0897    0.1429    0    0.1108]; 
elseif strcmp(campaign,'IMPACTS2022')
    %a = [1.107	1.370	0.98	1.037];
    %a_std = [0.046	0.083	0.014	0.128]; 
    % from 18 to 42°
    disp('Channel correction factors not final!')
    f_mirror = load('/Users/emma/Documents/PHIPS/PHIPS-HALO/Characterization/Pigtail/pigtail_mirror_fiber_coupling_test_correction_factors_220925.txt');
    f_pmt = load('/Users/emma/Documents/PHIPS/PHIPS-HALO/Characterization/MAPMT Crosstalk/pigtail_pmt_crosstalk_220923_gain75_angle.txt');
    f_pmt2 = load('/Users/emma/Documents/PHIPS/PHIPS-HALO/Characterization/MAPMT Crosstalk/pigtail_pmt_crosstalk_220926_gain75_angle.txt');
    f_pmt3 = load('/Users/emma/Documents/PHIPS/PHIPS-HALO/Characterization/MAPMT Crosstalk/pigtail_pmt_crosstalk_220928_gain75_angle.txt');
    f_pmt4 = load('/Users/emma/Documents/PHIPS/PHIPS-HALO/Characterization/MAPMT Crosstalk/pigtail_pmt_crosstalk_220929_gain75_angle.txt');
    Cm = f_mirror.*f_pmt;
    Cs = std([f_pmt; f_pmt2; f_pmt3; f_pmt4]);
elseif strcmp(campaign,'IMPACTS2023') % CHANGE THESE!
    % from 18 to 42°
    disp('Channel correction factors not defined!')
    Cm = [1	1	1	1];
    Cs = [0.029	0.029	0.054	0.021]; 
end


%% Read data
vbs = PhipsData.Properties.VariableNames;
Index_scattering = find(contains(vbs,'ScatteringAngle'));
Index_size = find(contains(vbs,'diameter'));
Index_area = find(contains(vbs,'proj_area'));
Index_ar = find(contains(vbs,'aspect_ratio'));
Index_T = find(contains(vbs,'Temperature'));
Index_RH = find(contains(vbs,'S_ice'));
Index_lat = find(contains(vbs,'lat'));
time_stamp=PhipsData.RealTimeStamp;
Start_time = PhipsData.RealTimeStamp(1);
end_time = PhipsData.RealTimeStamp(end);
ScaIn = table2array(PhipsData(:,...
    Index_scattering));
psize_c1c2 = table2array(PhipsData(:,Index_size));
parea_c1c2 = table2array(PhipsData(:,Index_area));
p_ar_c1c2 = table2array(PhipsData(:,Index_ar));
AR = min(p_ar_c1c2,[],2,'omitnan'); % smallest AR from the two cameras
parea = mean(parea_c1c2,2,'omitnan'); 
T = table2array(PhipsData(:,Index_T)); 
RH = table2array(PhipsData(:,Index_RH));
lat = table2array(PhipsData(:,Index_lat));
Tl = length(time_stamp);
sizeinfo=zeros(Tl,1);

% particle size larger than 26 micrometers from two cameras;
for ns = 1:Tl
    p1 = psize_c1c2(ns,1);
    p2 = psize_c1c2(ns,2);
    p1nan = isnan(p1);
    p2nan = isnan(p2);
    if p1nan && ~p2nan
        sizeinfo(ns)=p2;
    end
    if ~p1nan && p2nan
        sizeinfo(ns)=p1;
    end
    if p1nan && p2nan 
        sizeinfo(ns)=0;
    end
    if ~p1nan && ~p2nan
        sizeinfo(ns)=(p1+p2)/2;
    end
    if sizeinfo(ns)>0 && sizeinfo(ns)<sizelimit
        sizeinfo(ns)=0;
    end
end

% filtering by size;
ids = find(sizeinfo);
Sca_data=ScaIn(ids,:);
Time_data = time_stamp(ids,:);
Size_data = sizeinfo(ids);
T_data = T(ids);
RH_data = RH(ids);
AR_data = AR(ids);
lat_data = lat(ids);
area_data = parea(ids);

% filter NaN from the first two angles;
sca_temp = sum(Sca_data(:,1:2),2);
idds = find(~isnan(sca_temp));

% results;
time_s = Time_data(idds,:); % time
size_s = Size_data(idds,:); % size information
scain_s = Sca_data(idds,:); % scattering data
T_s = T_data(idds,:); % temperature
RH_s = RH_data(idds,:); % RH
lat_s = lat_data(idds,:); % latitude
AR_s = AR_data(idds,:); % aspect ratio
area_s = area_data(idds,:); % projected area
Nsignals = length(scain_s);

disp(['Altogether ',num2str(length(time_s)),' ice crystals were included into g analysis.'])


%% prepare the group average  
Ngroup = idivide(int64(Nsignals),int64(Ng));
remd = mod(Nsignals,Ng);
g_group = zeros(Ngroup,1);
C_group = g_group;
M_group = zeros(Ngroup,20);
T_mean_group = g_group;
T_std_group = g_group;
T_min_group = g_group;
T_max_group = g_group;
RH_mean_group = g_group;
RH_min_group = g_group;
RH_max_group = g_group;
RH_std_group = g_group;
lat_mean_group = g_group;
AR_mean_group = g_group;
AR_max_group = g_group;
AR_min_group = g_group;
AR_std_group = g_group;
time_group = time_s(1:Ngroup);
time_start_group = time_s(1:Ngroup);
time_end_group = time_s(1:Ngroup);
size_group = g_group;
area_group = g_group;
% get the group averaged intensity; 
for kg = 1:Ngroup 
    k1 = 1+(kg-1)*Ng;
    k2 = k1+Ng-1;
    % intensity; 
    sca_tmp = scain_s(k1:k2,:);
    M_group(kg,:) = nanmean(sca_tmp,1);
    % size 
    size_group(kg) = sum(size_s(k1:k2))/Ng;
    % area 
    area_group(kg) = sum(area_s(k1:k2))/Ng;
    % time
    dur = time_s(k2)-time_s(k1);
    time_group(kg) = time_s(k1)+dur/2;
    time_start_group(kg) = time_s(k1);
    time_end_group(kg) = time_s(k2);
    % Temperature
    T_mean_group(kg) = sum(T_s(k1:k2))/Ng;
    T_max_group(kg) = nanmax(T_s(k1:k2));
    T_min_group(kg) = nanmin(T_s(k1:k2));
    T_std_group(kg) = std(T_s(k1:k2),'omitnan');
    % RH
    RH_mean_group(kg) = sum(RH_s(k1:k2))/Ng;
    RH_max_group(kg) = nanmax(RH_s(k1:k2));
    RH_min_group(kg) = nanmin(RH_s(k1:k2));
    RH_std_group(kg) = std(RH_s(k1:k2),'omitnan');
    % Latitude
    lat_mean_group(kg) = sum(lat_s(k1:k2))/Ng;
    % AR
    AR_mean_group(kg) = sum(AR_s(k1:k2))/Ng;
    AR_max_group(kg) = nanmax(AR_s(k1:k2));
    AR_min_group(kg) = nanmin(AR_s(k1:k2));
    AR_std_group(kg) = std(AR_s(k1:k2),'omitnan');
end


%% ----- asymmetry factor ---- retrieval; 

% Correction factors according to normal distribution
% Cm is the mean and Cs the std correction factor
numt = 1000; % Number of samples 
cf = [];
for k = 1:length(Cm)
    cf(:,k) = Cm(k) + Cs(k).*randn(numt,1); 
end

% Asymmetry factor calculation
g = []; g_std = []; Cp = []; Cp_std = [];
for k = 1:size(M_group,1)
    temp = repmat(M_group(k,:),[numt 1]);
    gtt = []; Cptt = [];
    for kk = 1:length(Cm)
        temp(:,kk) = M_group(k,kk).*cf(:,kk);
    end
    [gtt, Cptt, ~] = asymmetryfactor_complexity(temp, repmat(size_group(k),[numt 1]));
    g(k) = mean(gtt,1,"omitnan");
    g_std(k) = nanstd(gtt,0,1);
    Cp(k) = mean(Cptt,1,"omitnan");
    Cp_std(k) = nanstd(Cptt,0,1);
end


%% --- Make data file ----
DATA_out = table(time_start_group,time_end_group, lat_mean_group, ...
    T_mean_group,T_max_group,T_min_group, ...
    RH_mean_group,RH_max_group,RH_min_group, ...
    AR_mean_group,AR_max_group,AR_min_group, ...
    size_group,g',g_std', Cp', Cp_std', 'VariableNames',...
    {'StartTime','EndTime', 'lat', ...
    'T_mean','T_max','T_min', ...
    'RHi_mean','RHi_max','RHi_min', ...
    'AR_mean','AR_max','AR_min', ...
    'diameter','g','g_std','Cp','Cp_std'});

% Calculate measured volume
delta_time = seconds(time_end_group-time_start_group);
delta_distance = 150.*delta_time;
delta_volume = 0.005./1e4.*delta_distance;
disp(['Average measurement volume ',num2str(nanmean(delta_volume)),' m3'])


end

