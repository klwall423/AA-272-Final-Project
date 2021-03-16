function [cosa_mean, cosa_w, edges] = calibration()
%% Read in relevant data

ephem = read_rinex_nav('abpo2000.10n');
%ephem = read_gps('ade12000.10n');
attitude = read_attitude('ICESAT_ATTITUDE_2010200.txt');
satpos = read_pos('ice12000.10pc');
gps = read_receiver('ice12000.10o');

%% Do Calibration

% find gps satellite positions
gpsstruct = struct('sat',[],'x',[],'y',[],'z',[]);
gpspos = struct('pos',gpsstruct);
for i = 1:length(gps.sec)
    % calculate measurement time in GPS seconds
    t = 86400 + gps.hour(i)*3600 + gps.min(i)*60 + gps.sec(i);
    sats = gps.gps(i).sat;
    gps_temp = struct('sat',[],'x',[],'y',[],'z',[]);
    for j = 1:length(sats)
        %pos = ephem2pos(ephem,sats(j),t);
        pos = ephem2pos2(ephem,sats(j),t);
        gps_temp.sat = [gps_temp.sat sats(j)];
        gps_temp.x = [gps_temp.x pos(1)];
        gps_temp.y = [gps_temp.y pos(2)];
        gps_temp.z = [gps_temp.z pos(3)];
    end
    gpspos.pos = [gpspos.pos gps_temp];
end
gpspos.pos(1) = [];

% calculate LOS vectors, SNRs, and attitudes
LOS = [];
SNR = [];
cosa = [];

for i=1:500%length(satpos.sec)
    % find corresponding gps data
    k = find(gps.hour == satpos.hour(i) & gps.min == satpos.min(i) & gps.sec == satpos.sec(i));
    if k
        gpspos_temp = [gpspos.pos(k).x',gpspos.pos(k).y',gpspos.pos(k).z'];
        satpos_temp = 10^3*[satpos.x(i)' satpos.y(i)' satpos.z(i)'];
        satpos_temp = repmat(satpos_temp,length(gpspos.pos(k).x),1);
        LOS_temp = bsxfun(@rdivide,(gpspos_temp-satpos_temp),vecnorm(gpspos_temp-satpos_temp,2,2));
        snr_temp = gps.gps(k).ss;
        SNR = [SNR snr_temp];
        [~,j] = min(abs(attitude.t_gps - (3600*satpos.hour(i)+60*satpos.min(i)+satpos.sec(i))));
        att_temp = [attitude.q1(j) attitude.q2(j) attitude.q3(j) attitude.q4(j)];
        
        % convert LOS to spacecraft frame
        dcm = dcmeci2ecef('IAU-2000/2006',[satpos.year(i) satpos.month(i) satpos.day(i) satpos.hour(i) satpos.min(i) satpos.sec(i)]);
        dcm2 = [cos(pi/4.2) sin(pi/4.2) 0;
               -sin(pi/4.2) cos(pi/4.2) 0;
               0 0 1];
        dcm = dcm2'*dcm';
        LOS_temp_eci = dcm*LOS_temp';   % LOS in ECI
        LOS_temp = quatrotate(att_temp,LOS_temp_eci');   % LOS in spacecraft frame
        LOS = [LOS; LOS_temp];
        ant = repmat([-1 0 0],length(LOS_temp),1);
        cosa_temp = dot(ant,LOS_temp,2);
        cosa = [cosa; cosa_temp];
        
    end
end

% sort SNR values into bins
[Nbins,edges,inds] = histcounts(SNR,100);
cosa_mean = zeros(length(Nbins),1);
cosa_w = zeros(length(Nbins),1);
for i = 1:length(Nbins)
    cosa_mean(i) = mean(cosa(inds==i));
    cosa_w(i) = std(cosa(inds==i));
end
mm = find(isnan(cosa_mean));
for nn = mm
    cosa_mean(nn) = cosa_mean(nn-1);
    cosa_w(nn) = cosa_w(nn-1);
end
pp = find(cosa_w == 0);
cosa_w(pp) = 0.0043;

figure;
scatter(SNR,cosa,'Filled');
xlabel('SNR');
ylabel('cos\alpha');


end