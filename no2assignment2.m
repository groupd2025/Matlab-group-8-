% Question 2: Group Member Information with Statistical Analysis
fprintf('=== COLLECTING AND ANALYZING GROUP MEMBER INFORMATION ===\n');

% Create a structure array to store group member information
group_members = struct();

% Add information for each group member (only first set of attributes)
% Member 1
group_members(1).name = 'John Doe';
group_members(1).background = 'Engineering';
group_members(1).home_district = 'Nairobi';
group_members(1).tribe = 'Kikuyu';
group_members(1).village = 'Gachie';
group_members(1).religion = 'Christian';
group_members(1).interests = {'Programming', 'Robotics', 'AI'};
group_members(1).age = 22;

% Member 2
group_members(2).name = 'Jane Smith';
group_members(2).background = 'Business';
group_members(2).home_district = 'Mombasa';
group_members(2).tribe = 'Luo';
group_members(2).village = 'Kisumu';
group_members(2).religion = 'Muslim';
group_members(2).interests = {'Marketing', 'Finance', 'Entrepreneurship'};
group_members(2).age = 23;

% Member 3
group_members(3).name = 'David Johnson';
group_members(3).background = 'Science';
group_members(3).home_district = 'Kakamega';
group_members(3).tribe = 'Luhya';
group_members(3).village = 'Webuye';
group_members(3).religion = 'Christian';
group_members(3).interests = {'Research', 'Mathematics', 'Astronomy'};
group_members(3).age = 24;

% Member 4
group_members(4).name = 'Sarah Wangari';
group_members(4).background = 'Arts';
group_members(4).home_district = 'Nakuru';
group_members(4).tribe = 'Kikuyu';
group_members(4).village = 'Naivasha';
group_members(4).religion = 'Christian';
group_members(4).interests = {'Writing', 'Poetry', 'Theater'};
group_members(4).age = 21;

% Convert to table for better visualization
member_table = struct2table(group_members);

% Display the table
disp('Group Members Information:');
disp(member_table);

% Save the variables to a MAT file
save('group_members_info.mat', 'group_members', 'member_table');

% Also save to Excel for easy viewing
writetable(member_table, 'group_members.xlsx');

fprintf('Group member information saved to group_members_info.mat and group_members.xlsx\n');

% Statistical Analysis and Visualization
fprintf('\n=== PERFORMING STATISTICAL ANALYSIS ===\n');

% Extract numerical data for analysis
ages = [group_members.age];

% Calculate basic statistics
age_stats = [mean(ages), median(ages), std(ages), min(ages), max(ages)];

% Display statistics
fprintf('Age Statistics (Mean, Median, Std, Min, Max): %.2f, %.2f, %.2f, %d, %d\n', age_stats);

% Create a summary table of statistics
stat_names = {'Mean', 'Median', 'Std', 'Min', 'Max'};
stats_table = table(age_stats', 'VariableNames', {'Age'}, 'RowNames', stat_names);
disp('Statistical Summary:');
disp(stats_table);

% Visualization 1: Demographic Information
figure('Position', [100, 100, 1200, 800]);

% Tribe distribution
subplot(2, 3, 1);
tribes = {group_members.tribe};
[unique_tribes, ~, tribe_idx] = unique(tribes);
tribe_counts = histcounts(tribe_idx, 1:length(unique_tribes)+1);
pie(tribe_counts);
title('Tribe Distribution');
legend(unique_tribes, 'Location', 'eastoutside');

% Religion distribution
subplot(2, 3, 2);
religions = {group_members.religion};
[unique_religions, ~, religion_idx] = unique(religions);
religion_counts = histcounts(religion_idx, 1:length(unique_religions)+1);
pie(religion_counts);
title('Religion Distribution');
legend(unique_religions, 'Location', 'eastoutside');

% Background distribution
subplot(2, 3, 3);
backgrounds = {group_members.background};
[unique_backgrounds, ~, background_idx] = unique(backgrounds);
background_counts = histcounts(background_idx, 1:length(unique_backgrounds)+1);
bar(background_counts);
set(gca, 'XTickLabel', unique_backgrounds);
title('Background Distribution');
ylabel('Number of Members');

% District distribution
subplot(2, 3, 4);
districts = {group_members.home_district};
[unique_districts, ~, district_idx] = unique(districts);
district_counts = histcounts(district_idx, 1:length(unique_districts)+1);
bar(district_counts);
set(gca, 'XTickLabel', unique_districts);
title('District Distribution');
ylabel('Number of Members');

% Age distribution
subplot(2, 3, 5);
histogram(ages, 'BinMethod', 'integers');
title('Age Distribution');
xlabel('Age');
ylabel('Number of Members');

% Village distribution
subplot(2, 3, 6);
villages = {group_members.village};
[unique_villages, ~, village_idx] = unique(villages);
village_counts = histcounts(village_idx, 1:length(unique_villages)+1);
bar(village_counts);
set(gca, 'XTickLabel', unique_villages);
title('Village Distribution');
ylabel('Number of Members');

sgtitle('Demographic Information of Group Members');
saveas(gcf, 'demographic_info.png');



% Visualization 3: Age vs Background
figure('Position', [100, 100, 800, 600]);

% Group ages by background
background_age_data = [];
background_labels = unique(backgrounds);
for i = 1:length(background_labels)
    background_ages = ages(strcmp(backgrounds, background_labels{i}));
    background_age_data = [background_age_data, background_ages'];
end


% Create a comprehensive report
fprintf('\n=== GENERATING COMPREHENSIVE REPORT ===\n');
report_filename = 'group_members_report.txt';
fid = fopen(report_filename, 'w');

fprintf(fid, 'COMPREHENSIVE GROUP MEMBERS REPORT\n');
fprintf(fid, '===================================\n\n');
fprintf(fid, 'Generated on: %s\n\n', datestr(now));

fprintf(fid, 'GROUP SUMMARY:\n');
fprintf(fid, 'Number of members: %d\n', length(group_members));
fprintf(fid, 'Average age: %.2f years\n', mean(ages));
fprintf(fid, 'Age range: %d to %d years\n\n', min(ages), max(ages));

fprintf(fid, 'DIVERSITY ANALYSIS:\n');
fprintf(fid, 'Tribes represented: %s\n', strjoin(unique_tribes, ', '));
fprintf(fid, 'Religions represented: %s\n', strjoin(unique_religions, ', '));
fprintf(fid, 'Backgrounds represented: %s\n', strjoin(unique_backgrounds, ', '));
fprintf(fid, 'Districts represented: %s\n', strjoin(unique_districts, ', '));
fprintf(fid, 'Villages represented: %s\n\n', strjoin(unique_villages, ', '));


fprintf(fid, 'DETAILED MEMBER PROFILES:\n');
fprintf(fid, '=========================\n\n');

for i = 1:length(group_members)
    fprintf(fid, 'MEMBER %d: %s\n', i, group_members(i).name);
    fprintf(fid, 'Background: %s\n', group_members(i).background);
    fprintf(fid, 'District: %s, Tribe: %s, Village: %s\n', ...
        group_members(i).home_district, group_members(i).tribe, group_members(i).village);
    fprintf(fid, 'Religion: %s\n', group_members(i).religion);
    fprintf(fid, 'Age: %d\n', group_members(i).age);
    % fprintf(fid, 'Interests: %s\n\n', strjoin(group_members(i).interests, ', '));
end

fclose(fid);
fprintf('Comprehensive report saved to: %s\n', report_filename);

% Create a summary table for presentation
summary_table = table(unique_tribes', histcounts(tribe_idx, 1:length(unique_tribes)+1)', ...
    'VariableNames', {'Tribe', 'Count'});
disp('Tribe Summary:');
disp(summary_table);

summary_table = table(unique_religions', histcounts(religion_idx, 1:length(unique_religions)+1)', ...
    'VariableNames', {'Religion', 'Count'});
disp('Religion Summary:');
disp(summary_table);

summary_table = table(unique_backgrounds', histcounts(background_idx, 1:length(unique_backgrounds)+1)', ...
    'VariableNames', {'Background', 'Count'});
disp('Background Summary:');
disp(summary_table);

fprintf('\n=== ALL TASKS COMPLETED ===\n');