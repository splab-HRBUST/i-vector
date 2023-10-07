fid = fopen('enroll.txt');
data = textscan(fid,'%s');
file_enroll = data{1,1};
save file_enroll file_enroll;