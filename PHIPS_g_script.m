function [DATA_out,scatt_data,legcoefs_m] = PHIPS_g_script(PhipsData, Ng, sizelimit, campaign)
% Calculates g for a group of particles
% Output:   DATA_out
%           scatt_data ->  Average scattering function used for g retrieval
%           legcoefs_m ->  Legendre coefficients
% Input:    PhipsData  ->  PHIPS level 5 table, data not corrected
%           Ng         ->  Number of particles per group
%           sizelimit  ->  Lower size limit, e.g. 15 micron
%           campaign   ->  Name of the campaign
% Note:     set Ng = [] if using one averaged ASF
addpath('/Users/emma/Documents/m-codes/PHIPS analysis/AsymmetryFactorAnalysis/src/chebfun-master');


%% Channel correction factor std
if strcmp(campaign,'ACLOUD') || strcmp(campaign,'SOCRATES')
    % Based on individual droplets measured during ARISTO2017
    % from 18 to 42°
    Cm = [1.0313    0.9715    1.3200    1.0924    2.0480    1.3183    1.8443];
    Cs = [0.2656    0.1490    0.2057    0.1718    0.4587    0.3283    0.4339]; 
elseif strcmp(campaign,'CIRRUS-HL') 
    % Glass bead calibration results performed during the campaign 
    % 20 micron beads
    % Uses residual method
    % Laser polarisation measured 50-50
    % updated on 28.06.2024
    % from 18 to 90°
    %Cm = [1.0434    0.9651    1.0858    0.8650    1.0681    1.1828    0.8938    1.2364    0.9245    0.9800];
    %Cs = [0.0664    0.0342    0.0740    0.0741    0.0989    0.1339    0.0881    0.2565    0.1004    0.1598]; 
    % from 18 to 66° (old, laser polarisation 25-75)
    %Cm = [1.02770721259156  0.950933977663321   1.09191594324931    0.810532702269994   1.11678247431606    1.14459909547292    1.01252505501933    1   1   1   1   1   1   1   1   1   1   1   1   1];
    %Cs = [0.046682321334456 0.0393249176506893  0.0729291837694767  0.0701617985969     0.112478812681935   0.144840032467046   0.10673100846613    0.286773108394354   0.112565305479333   0.166970805882387   0.119614380630246   0.127862865786036   0.1281710441128     0.102649522076128   0.158685538675458   0.134215587987003   0.125794543516891   0.150365961683902   0.126177997719152   0.241783006848424]; 
    % from 18 to 66°
    Cm = [1.03262044205187  0.945182324295287   1.08590099259698    0.77668419817837    1.03477253590165    1.07226582120809    0.905489412301888   1   1   1   1   1   1   1   1   1   1   1   1   1];
    Cs = [0.0471821174028975,0.039095831291556,0.072259419508463,0.0670821161859991,0.104336240387173,0.135199107712607,0.0949887576812448,0.287450922882544,0.104275678253414,0.158365102736563,0.120433820193251,0.144922071511689,0.151589833635204,0.117860141660025,0.175556018058955,0.150254043287574,0.13024568930038,0.139583452864309,0.143987169410693,0.318276605714386]; 
elseif strcmp(campaign,'IMPACTS2022')
    % Glass bead calibration results performed during the campaign
    % from 18 to 82°
    % updated on 12.12.2023
    Cm = [1.0363    1.2316    0.8692    1.0445     0.9543   0.9249  0.9825  1.1466	 0.8989];
    Cs = [0.0966    0.0775    0.0571    0.1470     0.0834   0.1011  0.0734  0.2365   0.1772]; 
elseif strcmp(campaign,'IMPACTS2023') % CHANGE THESE!
    % from 18 to 42°
    disp('Channel correction factors not defined!')
    Cm = [1	1	1	1];
    Cs = [0.029	0.029	0.054	0.021]; 
elseif strcmp(campaign,'RICE04')
    % Glass bead calibration results performed before the campaign
    % from 6° to 166°
    % updated on 24.04.2024
    Cm = [2.1295 1.1180 1.0335 0.99084 0.96545 1.0678 1.0438 1.0144 0.7163 1.0501 0.8210 0.8936 0.7488 0.9496 1.1200 1.0230 1.0610 0.9967 1.0909 0.8386 1.5101];
    Cs = [0.3386 0.1134 0.1140 0.0835 0.1030 0.1294 0.1798 0.1982 0.1218 0.4210 0.1727 0.1952 0.1576 0.2310 0.3060 0.5161 0.4855 0.3875 0.3948 0.2204 0.2886];
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
if ~isempty(Index_T) % if level 5 data exists
    T = table2array(PhipsData(:,Index_T)); 
    RH = table2array(PhipsData(:,Index_RH));
    lat = table2array(PhipsData(:,Index_lat));
else
    T = NaN.*psize_c1c2; 
    RH = NaN.*psize_c1c2;
    lat = NaN.*psize_c1c2;
end
Tl = length(time_stamp);
sizeinfo = zeros(Tl,1);

% particle size larger than sizelimit from two cameras;
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
Sca_data = ScaIn(ids,:);
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
if isempty(Ng) % only one g,Cp retrieval
    Ngroup = 1;
    Ng = Nsignals;
else
    Ngroup = idivide(int64(Nsignals),int64(Ng));
end
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
g = []; g_std = []; Cp = []; Cp_std = []; scatt_data = []; C0 = []; C0_std = [];
for k = 1:size(M_group,1)
    temp = repmat(M_group(k,:),[numt 1]);
    gtt = []; Cptt = [];
    for kk = 1:length(Cm)
        temp(:,kk) = M_group(k,kk).*cf(:,kk);
    end
    [gtt, Cptt, C0tt, legcoefs_m] = asymmetryfactor_complexity(temp, repmat(size_group(k),[numt 1]));
    g(k) = mean(gtt,1,"omitnan");
    g_std(k) = nanstd(gtt,0,1);
    Cp(k) = mean(Cptt,1,"omitnan");
    Cp_std(k) = nanstd(Cptt,0,1);
    C0(k) = mean(C0tt,1,"omitnan");
    C0_std(k) = nanstd(C0tt,0,1);
    scatt_data(k,:) = mean(temp,1,"omitnan");
end


%% --- Make data file ----
DATA_out = table(time_start_group,time_end_group, lat_mean_group, ...
    T_mean_group,T_max_group,T_min_group, ...
    RH_mean_group,RH_max_group,RH_min_group, ...
    AR_mean_group,AR_max_group,AR_min_group, ...
    size_group,g',g_std', Cp', Cp_std', C0', C0_std', ...
    'VariableNames',...
    {'StartTime','EndTime', 'lat', ...
    'T_mean','T_max','T_min', ...
    'RHi_mean','RHi_max','RHi_min', ...
    'AR_mean','AR_max','AR_min', ...
    'diameter','g','g_std','Cp','Cp_std', ...
    'C0','C0_std'});

% Calculate measured volume
delta_time = seconds(time_end_group-time_start_group);
delta_distance = 150.*delta_time;
delta_volume = 0.005./1e4.*delta_distance;
disp(['Average measurement volume ',num2str(nanmean(delta_volume)),' m3'])


end

