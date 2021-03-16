clc
clear
format long g

GM = 3.986005e14;             % earth's universal gravitational [m^3/s^2]
c = 2.99792458e8;             % speed of light (m/s)
omegae_dot = 7.2921151467e-5; % earth's rotation rate (rad/sec)

% read navigation file to get orbit elements
%ephemeris = read_rinex_nav('brdc2000.06n');
ephemeris = read_rinex_nav('abpo2000.10n');

tsv = 86400; % space vehicle time [s]

% initialize constants and variables
svid = ephemeris(:,1);
m0   = ephemeris(:,2);
dn   = ephemeris(:,3);
e    = ephemeris(:,4);
a    = (ephemeris(:,5)).^2;
omg0 = ephemeris(:,6);
i0   = ephemeris(:,7);
w    = ephemeris(:,8);
odot = ephemeris(:,9);
idot = ephemeris(:,10);
cuc  = ephemeris(:,11);
cus  = ephemeris(:,12);
crc  = ephemeris(:,13);
crs  = ephemeris(:,14);
cic  = ephemeris(:,15);
cis  = ephemeris(:,16);
toe  = ephemeris(:,17);
iode = ephemeris(:,18);
GPS_week = ephemeris(:,19);
toc=ephemeris(:,20);
af0= ephemeris(:,21);
af1= ephemeris(:,22);
af2= ephemeris(:,23);
TGD=ephemeris(:,24);

nn = length(ephemeris);
for ii = 1:nn
    % Procedure for coordinate calculation
    n0 = sqrt(GM/a(ii)^3); % (rad/s)    
    tk = tsv-toe(ii);      % Time from eph ref epoch (s)
    n = n0+dn(ii);         % Corrected mean motion (rad/s)
    M = m0(ii)+n*tk;       % Mean anomaly (rad/s)
    
    % Perform Newton-Raphson solution for Eccentric anomaly estimate (rad)
	NRnext = 0;
	NR = 1;
	m = 1;
    while abs(NRnext-NR)>1e-15;
        NR=NRnext;
        f=NR-e(ii)*sin(NR)-M;
        f1=1-e(ii)*cos(NR);
        f2=e(ii)*sin(NR);
        NRnext=NR-(f/(f1-(f2*f/2*f1)));
        m=m+1;
    end
	E=NRnext; % Eccentric anomaly estimate for computing delta_tr (rad)
    
    % Time correction
    F = -2*sqrt(GM)/c^2; % (s/m^1/2)
    delta_tr = F*e(ii)*sqrt(a(ii))*sin(E);
    delta_tsv = af0(ii)+af1(ii)*(tsv-toe(ii))+delta_tr;
    t = tsv-delta_tsv;
    
    tk=t-toe(ii);  		 % Time from eph ref epoch (s)
    M=m0(ii)+n*tk;	     % Mean anomaly (rad/s)
    
    % Perform Newton-Raphson solution for Eccentric anomaly (rad)
    NRnext=0;
    NR=1;
    m=1;
    while abs(NRnext-NR)>1e-15;
        NR=NRnext;
        f=NR-e(ii)*sin(NR)-M;
        f1=1-e(ii)*cos(NR);
        f2=e(ii)*sin(NR);
        NRnext=NR-(f/(f1-(f2*f/2*f1)));
        m=m+1;
    end;
    E=NRnext; % Eccentric anomaly (rad)

    v = atan2(sqrt(1-e(ii)^2)*sin(E), cos(E)-e(ii));
    phi = v+w(ii);
    u = phi                    + cuc(ii)*cos(2*phi)+cus(ii)*sin(2*phi);
    r = a(ii)*(1-e(ii)*cos(E)) + crc(ii)*cos(2*phi)+crs(ii)*sin(2*phi);
    i = i0(ii)+idot(ii)*tk     + cic(ii)*cos(2*phi)+cis(ii)*sin(2*phi);
    Omega = omg0(ii)+(odot(ii)-omegae_dot)*tk-omegae_dot*toe(ii);
    x1 = cos(u)*r;
    y1 = sin(u)*r;
	
	% ECEF coordinates
    satp(1,ii) = svid(ii);
    satp(2,ii) = x1*cos(Omega)-y1*cos(i)*sin(Omega);
    satp(3,ii) = x1*sin(Omega)+y1*cos(i)*cos(Omega);
    satp(4,ii) = y1*sin(i);
end

