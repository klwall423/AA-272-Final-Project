function rec_data = read_receiver(rec_file)

% initialize data structure
gpsstruct = struct('sat',[],'pr',[],'ss',[]);
rec_data = struct('year',[],'month',[],'day',[],'hour',[],'min',[],'sec',[],'nsats',[],'gps',gpsstruct);

file = fopen(rec_file);
% skip header lines
for i=1:15
    fgetl(file);
end

i=1;
while ~feof(file)
    % collect time data
    timedata = fgetl(file);
    rec_data.year = [rec_data.year str2num(timedata(2:3))];
    rec_data.month = [rec_data.month str2num(timedata(5:6))];
    rec_data.day = [rec_data.day str2num(timedata(8:9))];
    rec_data.hour = [rec_data.hour str2num(timedata(11:12))];
    rec_data.min = [rec_data.min str2num(timedata(14:15))]; %#ok<*ST2NM>
    rec_data.sec = [rec_data.sec str2double(timedata(17:26))];
    
    nsat = str2num(timedata(32));    % get number of satellites in time period
    rec_data.nsats = [rec_data.nsats nsat];
    
    gps_temp = struct('sat',[],'pr',[],'ss',[]);
    
    % collect satellite numbers
    jj = 34;    % initialize
    for kk = 1:nsat
        gps_temp.sat = [gps_temp.sat str2num(timedata(jj:jj+1))];
        jj = jj+3;
    end
    
    % collect receiver data
    for kk = 1:nsat
        gpsdata = fgetl(file);
        gps_temp.pr = [gps_temp.pr str2double(gpsdata(35:48))];
        gpsdata = fgetl(file);
        gps_temp.ss = [gps_temp.ss str2double(gpsdata(23:32))];
    end
rec_data.gps = [rec_data.gps gps_temp];
i=i+1;
end

% remove extra gps entry
rec_data.gps(1) = [];

end