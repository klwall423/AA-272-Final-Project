function att_data = read_attitude(att_file)

% initialize structure
att_data = struct('t_gps',[],'q1',[],'q2',[],'q3',[],'q4',[]);

% read data
nums = readmatrix(att_file);
% assign data to structure
att_data.t_gps = nums(:,1);
att_data.q1 = nums(:,2);
att_data.q2 = nums(:,3);
att_data.q3 = nums(:,4);
att_data.q4 = nums(:,5);

end