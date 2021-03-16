clear;
%% Read in relevant data

ephem = read_rinex_nav('abpo2000.10n');
attitude = read_attitude('ICESAT_ATTITUDE_2010200.txt');
satpos = read_pos('ice12000.10pc');
gps = read_receiver('ice12000.10o');

%% Extract data

[cosa_mean, cosa_w, edges] = calibration();

% find gps satellite positions
gpsstruct = struct('sat',[],'x',[],'y',[],'z',[]);
gpspos = struct('pos',gpsstruct);
for i = 1:length(gps.sec)
    % calculate measurement time in GPS seconds
    t = (20-gps.day(i))*86400 + gps.hour(i)*3600 + gps.min(i)*60 + gps.sec(i);
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
LOS_ECI = [];
SNR = [];
att = [];
times = [];

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
        [bb,j] = min(abs(attitude.t_gps - (3600*satpos.hour(i)+60*satpos.min(i)+satpos.sec(i))));
        att_temp = [attitude.q1(j) attitude.q2(j) attitude.q3(j) attitude.q4(j)];
        
        % convert LOS to spacecraft frame
        dcm = dcmeci2ecef('IAU-2000/2006',[satpos.year(i) satpos.month(i) satpos.day(i) satpos.hour(i) satpos.min(i) satpos.sec(i)]);
        dcm2 = [cos(pi/4.2) sin(pi/4.2) 0;
               -sin(pi/4.2) cos(pi/4.2) 0;
               0 0 1];
        dcm = dcm2'*dcm';

        LOS_temp_eci = dcm*LOS_temp';   % LOS in ECI
        LOS_ECI = [LOS_ECI; LOS_temp_eci'];
        satpos_eci_temp = dcm*satpos_temp';

        LOS_temp = quatrotate(att_temp,LOS_temp_eci');   % LOS in spacecraft frame
        LOS = [LOS; LOS_temp];
        att_temp = repmat(att_temp,[length(LOS_temp),1]);
        att = [att; att_temp];
        time_temp = (gps.hour(k)*3600 + gps.min(k)*60 + gps.sec(k))*ones(1,length(LOS_temp));
        times = [times time_temp];
        ant = repmat([-1 0 0],length(LOS_temp),1);
    end
end

%% Attitude Estimation

% use first N measurements
N = 500;
bls = zeros(N,3);
x_ref = zeros(N,3);
vs = zeros(N,3);
vs_nf = zeros(N,3);
qs = zeros(N,4);
q_ref = zeros(N,4);

count = 1;
for p = 1:N
    k = find(times == times(count));
    l = length(k);
    % simulate additional vector measurements
    v_sensor = [0 sqrt(2)/2 sqrt(2)/2];
    vs(p,:) = quatrotate(quatinv(att(k(1),:)),v_sensor)+normrnd(0,0.01,1,3);
    vs_nf(p,:) = quatrotate(quatinv(att(k(1),:)),v_sensor);
    
    x_ref(p,:) = quatrotate(quatinv(att(k(1),:)),[-1 0 0]);
    q_ref(p,:) = att(k(1),:);

    c = zeros(l,1);
    W = zeros(l,l);
    H = zeros(l,3);
    for i=1:l
        % look up SNR to find corresponding cosa
        j = find(SNR(k(i)) < edges);
        if isempty(j)
            j = length(edges);
        end
        c(i) = cosa_mean(j(1)-1);
        % weights matrix
        W(i,i) = 1/(cosa_w(j(1)-1))^2;
        % H matrix
        H(i,:) = LOS_ECI(k(i),:);
    end

    % estimate bl0
    b0l = inv(H'*W*H)*H'*W*c;
    b0l = b0l/norm(b0l);

    % Solve for bl
    theta = acos(b0l(3));
    phi = acos(b0l(1)/sin(theta));
    x0 = [theta; phi];%b0l;
    J = @(x)1/2*sum(W*abs(c-H*[sin(x(1))*cos(x(2)); sin(x(1))*sin(x(2)); cos(x(1))]).^2);
    x = fminsearch(J,x0);
    bls(p,:) = [sin(x(1))*cos(x(2)) sin(x(1))*sin(x(2)) cos(x(1))];
    %bls(p,:) = x';
    count = k(end)+1;
    %bls(p,:) = b0l;
    
    % implement TRIAD method
    A = triad([-1 0 0]',[0 sqrt(2)/2 sqrt(2)/2]',bls(p,:)',vs(p,:)');  
    Ap = triad([-1 0 0]',[0 sqrt(2)/2 sqrt(2)/2]',x_ref(p,:)',vs_nf(p,:)');
    qs(p,:) = dcm2quat(A);
    qps(p,:) = dcm2quat(Ap);
    %count = count + 1;
end

eul_ref = rad2deg(quat2eul(qps));
eul = rad2deg(quat2eul(qs));

%k = find(eul(:,1) > 0);
%eul(k,1) = eul(k,1) - 360;

figure;
plot(eul_ref(:,1),'LineWidth',3);
hold on;
scatter(linspace(1,N,N),eul(:,1),'Filled');
xlabel('Time (s)');
ylabel('\psi (deg)');
legend('Actual','Estimate');

figure;
plot(eul_ref(:,2),'LineWidth',3);
hold on;
scatter(linspace(1,N,N),eul(:,2),'Filled');
xlabel('Time (s)');
ylabel('\theta (deg)');
legend('Actual','Estimate');

figure;
plot(eul_ref(:,3),'LineWidth',3);
hold on;
scatter(linspace(1,N,N),eul(:,3),'Filled');
xlabel('Time (s)');
ylabel('\phi (deg)');
legend('Actual','Estimate');

err = eul - eul_ref;
k = find(err < -350);
err(k) = err(k) + 360;

figure;
scatter(linspace(1,N,N),err(:,1),'Filled');
hold on;
scatter(linspace(1,N,N),err(:,2),'Filled');
scatter(linspace(1,N,N),err(:,3),'Filled');
xlabel('Time (s)');
ylabel('Error (deg');
legend('\Delta \psi','\Delta \theta','\Delta \phi');

figure;
plot3(x_ref(:,1),x_ref(:,2),x_ref(:,3),'LineWidth',3);
hold on;
scatter3(bls(:,1),bls(:,2),bls(:,3),'Filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
ylim([-1 1]);

vecerr = bls-x_ref;
figure;
scatter(linspace(1,N,N),vecerr(:,1),'Filled');
hold on;
scatter(linspace(1,N,N),vecerr(:,2),'Filled');
scatter(linspace(1,N,N),vecerr(:,3),'Filled');
xlabel('Time (s)');
ylabel('Error');
legend('\Delta X','\Delta Y','\Delta Z');
