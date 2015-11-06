%GUINEA

% import .csv value
filename= 'matlab_data.csv';
%%%%%%%%%%%%%%% MODIFY FILE PATH TO DATA DESIRED ABOVE %%%%%%%%%%%%%%% 
fid = fopen(filename);
raw_data=textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
fclose(fid);

% data of the form [day of outbreak, cases (i.e. total infected) on this day, deaths
% (i.e. removed) by this day]
lib_data2 = [raw_data{2}, raw_data{3}, raw_data{11}];
temp_matrix = zeros(size(lib_data2));

count=1;
for i=1:length(lib_data2)
    lib_data2(i,1) = lib_data2(i,1); %- 102; % shift the data to start on this
                                           % (arbitrary!) date
    if (~isnan(lib_data2(i,2)) & ~isnan(lib_data2(i,3)))
        if (lib_data2(i,1) >= 0)
            temp_matrix(count,:) = lib_data2(i,:);
            count=count+1;
        end
    end
    
end
lib_data=flipud(temp_matrix(find(temp_matrix(:,1),1,'first'):find(temp_matrix(:,1),1,'last')+1,:));

%Adjusting total cases to infected only
lib_data(:,2) = lib_data(:,2) - lib_data(:,3);

A = lib_data(2:72,:);
B = lib_data(74:end,:);
lib_data=[A;B]

%Adjusting from cumulative case counts to current
% for i = 2:length(lib_data(:,2))
%     temp(i-1) = max(0,lib_data(i,2) - lib_data(i-1,2));
% end
% lib_data(2:end,2) = temp;

tSpan = lib_data(:,1);

plot(lib_data(:,1),lib_data(:,2),'b.',lib_data(:,1), lib_data(:,3), 'rs','MarkerSize',10)
legend('Cumulative non-death cases data','Deaths data','Location','Northwest')

