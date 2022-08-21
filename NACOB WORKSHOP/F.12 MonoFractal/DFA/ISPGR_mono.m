%% TITLE: ISPGR_mono.m
% DATE: June, 2022
% AUTHOR: Tyler M. Wiles, MS
% EMAIL: tylerwiles@unomaha.edu

% dfa.m function questions should be directed toward bmchnonan@unomaha.edu

% UPDATED: June 2022

% DESCRIPTION:
% This m-file accompanies the ISPGR 2022 workshop titled "Univeriate and
% bivariate fractal methods for movement science." Additionally, this file
% works in tandem with Part 1 of the lecture slides on Detrended
% Fluctuation Analysis.
    % Section 1 - DFA analysis on one walking timeseries
    % Section 2 - DFA analysis on multiple walking timeseries
    % Section 3 - DFA analysis on one COP timeseries
    % Section 4 - DFA analysis on multiple COP timeseries

% Redistribution and use of this script, with or without
% modification, is permitted provided this copyright notice,
% the original authors name and the following disclaimer remains.

% DISCLAIMER: It is the user's responsibility to check the code is returning
% the appropriate results before reporting any scientific findings.
% The author will not accept any responsibility for incorrect results
% directly relating to user error.

% Required Toolboxes:
% Statistics and Machine Learning Toolbox
% Signal Processing Toolbox
% Image Processing Toolbox

%% Section 1: DFA analysis on one stride interval timeseries

% This section of the tutorial will briefly show how to run DFA on stride
% interval data from a single participant while walking at their self-selected
% pace overground. 

% Clearing everything to start from scratch
clear all; close all; clc;

% This line will allow us to use the dfa.m function in the FUNCTION folder
addpath('FUNCTION')

% Loading in data. Column 1 is time, Column 2 is stride intervals.
s206_selfpaced = readmatrix('s206_selfpaced.csv');

%% Defining our input parameters for the dfa.m function

% Select column 2 from our data which is stride interval
ts = s206_selfpaced(:,2);

% n_min is set to 16 as recommended by papers in the lecture notes of Part 1
% slide 27. Citations can be found on line 26 of the dfa.m function.
n_min = 16;

% n_max set to a default of the length of the time series divided by 9
% as recommended by Damouras, et al., 2010. The appropriate paper can be
% found on line 25 of the dfa.m function. This input is also discussed
% in the lecture on Part 1 Slide 27.
n_max = length(ts)/9;

% n_length is set to the default number of points to sample best fit.
n_length = 18;

% The next line of code runs the actual DFA analysis. This function runs
% through the steps outlined in the lecture, specifically Part 1 slides
% 18-23. The output, a, represents the alpha value and is what we are most
% interested in when comparing between subjects, trials, conditions etc.
% Take note of how the alpha value compares to the DFA analysis
% interpretation slide in Part 1 slides 24-25. Is it persistent or
% antipersistent? You will also produce a figure that can be interpreted
% using the slides 24-25 from Part 1 of the lecture.

[a, r2] = dfa(ts, n_min, n_max, n_length)
title('DFA Stride Intervals s206');
xlabel('log(n)');
ylabel('log(F(n))');
set(gca, 'FontSize', 18); % Setting Font Size
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

% Done!

%% Section 2: DFA analysis on multiple stride interval timeseries

% Clearing everything to start from scratch
clear all; close all; clc;

% This line will allow us to use the dfa.m function in the FUNCTION folder
addpath('FUNCTION')

% Setup for the loop

% Gets directory to run the rest of the code. Simply choose where you want
% the data to be stored and run.
my_directory = uigetdir;

% This is the folder where our results will go
output_directory = append(my_directory, '/ANALYSIS OUTPUT');

% Identify all files with selfpaced.csv, which is our stride interval data
my_files = dir(fullfile(my_directory, '*selfpaced.csv'));

dfa_values = []; %This is where our DFA reults will go

for k = 1:length(my_files)
    base_filename = my_files(k).name;
    full_filename = fullfile(my_directory, base_filename);
    
    % Loading in data. Column 1 is time, Column 2 is stride intervals.
    temp_dat = readmatrix(full_filename);
    
    %% Defining our input parameters for the dfa.m function
    
    % Select column 2 from our data corresponding to stride interval
    ts = temp_dat(:,2);
    
    % n_min is set to 16 as recommended by papers in the lecture notes of Part 1
    % slide 27. Citations can be found on line 26 of the dfa.m function.
    n_min = 16;
    
    % n_max set to a default of the length of the time series divided by 9
    % as recommended by Damouras, et al., 2010. The appropriate paper can be
    % found on line 25 of the dfa.m function. This input is also discussed
    % in the lecture on Part 1 Slide 27.
    n_max = length(ts)/9;
    
    % n_length is set to the default number of points to sample best fit.
    % This input can be moved outside the loop if preferred.
    n_length = 18;
    
    % Up until this point the code is not much different than Section 1. Now we
    % will produce a single figure with our time series plotted, our DFA
    % output, and a histogram of our time series.
    
    %% Plot and run DFA
    
    % Figure 1: Stride Intervals overtime
    subplot(2, 2, [1,2]);
    plot(ts);
    title(base_filename);
    xlabel('Time (s)');
    ylabel('Stride Interval (seconds)');
    set(gca, 'FontSize', 15); % Setting Font Size
    
    % Figure 2: DFA (Exact same process as Section 1)
    subplot(2, 2, 3);
    [a, r2] = dfa(ts, n_min, n_max, n_length)
    title('DFA Stride Intervals');
    xlabel('log(n)');
    ylabel('log(F(n))');
    set(gca, 'FontSize', 15); % Setting Font Size
    
    % Figure 3: Histogram of stride intervals overtime
    subplot(2, 2, 4);
    histogram(ts);
    title('Stride Interval Distribution');
    xlabel('Stride Interval (seconds)');
    ylabel('Frequency of Occurance');
    set(gca, 'FontSize', 15); % Setting Font Size
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
    % The next two lines of code will require you to press any key on the
    % keyboard to close the figure and progress to the next trial
    disp('Press any key to continue!')
    pause;
    
    %% Data collating for analysis
    
    % Extracting filenames for exporting
    filename_split = strsplit(base_filename, '.');
    png_filename = filename_split{1};
    
    % Collating our figures into the folder titled "ANALYSIS OUTPUT"
    png = append(png_filename, '.png');
    saveas(gcf,fullfile(output_directory, png));
    
    % Collating alpha values from our DFA function with corresponding
    % subject id and condition.
    dfa_values{k,1} = png_filename;
    dfa_values{k,2} = a;
    
    close all % Close figures before looping through to the next trial
end

% Exporting DFA values into a table in long format for statistical
% analysis in R.
dfa_results = cell2table(dfa_values, 'VariableNames', {'id','alpha'});
writetable(dfa_results, fullfile(output_directory, 'Stride Interval DFA Values.csv'));

% Done!

%% Section 3: DFA analysis on one COP timeseries

% Center of Pressure data for Subject 28 of the Human Balance Evaluation
% Database (https://archive.physionet.org/physiobank/database/hbedb/)
% The subject performed a 60 second standing task with their eyes closed
% on a firm surface.

% This section of the tutorial will briefly show how to run DFA on center
% of pressure velocity for a single subject for a single trial. Notice,
% this procedure is the exact same as Section 1.

% Clearing everything to start from scratch
clear all; close all; clc;

% This line will allow us to use the dfa.m function in the FUNCTION folder
addpath('FUNCTION')

% Loading in data. Column 1 is time, Column 3 is COPy Velocity.
s018_cop = readmatrix('s018_cop.csv');

%% Defining our input parameters for the dfa.m function

% Select column 3 from our data which is COPy velocity
ts = s018_cop(:,3);

% n_min is set to 16 as recommended by papers in the lecture notes of Part 1
% slide 27. Citations can be found on line 26 of the dfa.m function.
n_min = 16;

% n_max set to a default of the length of the time series divided by 9
% as recommended by Damouras, et al., 2010. The appropriate paper can be
% found on line 25 of the dfa.m function. This input is also discussed
% in the lecture on Part 1 Slide 27.
n_max = length(ts)/9;

% n_length is set to the default number of points to sample best fit.
n_length = 18;

% The next line of code runs the actual DFA analysis. This function runs
% through the steps outlined in the lecture, specifically Part 1 slides
% 18-23. The output, a, represents the alpha value and is what we are most
% interested in when comparing between subjects, trials, conditions etc.
% Take note of how the alpha value compares to the DFA analysis
% interpretation slide in Part 1 slides 24-25. Is it persistent or
% antipersistent? You will also produce a figure that can be interpreted
% using the slides 24-25 from Part 1 of the lecture.

[a, r2] = dfa(ts, n_min ,n_max, n_length)
title('DFA COPy Velocity s018');
xlabel('log(n)');
ylabel('log(F(n))');
set(gca, 'FontSize', 18); % Setting Font Size
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

% Dealing with the crossover
[a, r2, out_a, out_l] = dfa(ts, n_min ,n_max, n_length);
title('DFA COPy Velocity s018');
xlabel('log(n)');
ylabel('log(F(n))');
set(gca, 'FontSize', 18); % Setting Font Size
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

% Adjusting regression lines by choosing our breakpoint
breakpoint= 98;

n = out_a(:,1); % box sizes
F = out_a(:,2); % fluctuation for a given box size
logn = log(n); % log box sizes
logF = log(F); % log fluctuation for a given box size

fontSize = 15;
markerSize = 20;     
[ashort,along] = crossover(logn, logF, breakpoint, 1)
title('DFA COPy Velocity short and long fluctuation');
xlabel('log(n)');
ylabel('log(F(n))');
caption = sprintf('long = %.4f , short = %.4f ', along, ashort);
text(min(logn), max(logF),caption, 'FontSize', fontSize);

% Done!

%% Section 4: DFA analysis on multiple COP timeseries

% Clearing everything to start from scratch
clear all; close all; clc;

% Setup for the loop

% Gets directory to run the rest of the code. Simply choose where you want
% the data to be stored and run.
my_directory = uigetdir;

% This is the folder where our results will go
output_directory = append(my_directory, '/ANALYSIS OUTPUT');

% Identify all files with cop.csv, which is our center of pressure velocity data
my_files = dir(fullfile(my_directory, '*cop.csv'));

dfa_values = []; % This is where our DFA reults will go

for k = 1:length(my_files)
    base_filename = my_files(k).name;
    full_filename = fullfile(my_directory, base_filename);
    
    % Loading in data. Column 1 is time, Column 3 is COPy Velocity
    temp_dat = readmatrix(full_filename);
    
    %% Defining our input parameteres for the dfa.m function
    
    % Select column 2 from our data corresponding to velocity
    ts = temp_dat(:,3);
    
    % n_min is set to 16 as recommended by papers in the lecture notes of Part 1
    % slide 27. Citations can be found on line 26 of the dfa.m function.
    n_min = 16;
    
    % n_max set to a default of the length of the time series divided by 9
    % as recommended by Damouras, et al., 2010. The appropriate paper can be
    % found on line 25 of the dfa.m function. This input is also discussed
    % in the lecture on Part 1 Slide 27.
    n_max = length(ts)/9;
    
    % n_length is set to the default number of points to sample best fit.
    % This input can be moved outside the loop if preferred.
    n_length = 18;
    
    % Up until this point the code is not much different than Section 1. Now we
    % will produce a single figure with our time series plotted, our DFA
    % output, and a histogram of our time series.
    
    %% Plot and run DFA
    
    % Figure 1: COPy velocity overtime
    subplot(2, 2, [1,2]);
    plot(ts);
    title(base_filename);
    xlabel('Time (s)');
    ylabel('Velocity (cm/s)');
    set(gca, 'FontSize', 15); %Setting Font Size
    
    % Figure 1: DFA (Exact same process as Section 1)
    subplot(2, 2, 3);
    [a, r2] = dfa(ts, n_min, n_max, n_length)
    title('DFA COPy Velocity');
    xlabel('log(n)');
    ylabel('log(F(n))');
    set(gca, 'FontSize', 15); % Setting Font Size
    
    % Figure 3: Histogram of COPy velocity overtime
    subplot(2, 2, 4);
    histogram(ts);
    title('Velocity Distribution');
    xlabel('Velocity (cm/s)');
    ylabel('Frequency of Occurance');
    set(gca, 'FontSize', 15); % Setting Font Size
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
    % The next two lines of code will require you to press any key on the
    % keyboard to close the figure and progress to the next trial
    disp('Press any key to continue!')
    pause;
    
    %% Data collating for analysis
    
    % Extracting filenames for exporting
    filename_split = strsplit(base_filename, '.');
    png_filename = filename_split{1};
    
    % Collating our figures into the folder titled "ANALYSIS OUTPUT"
    png = append(png_filename, '.png');
    saveas(gcf,fullfile(output_directory, png));
    
    % Collating alpha values from our DFA function with corresponding
    % subject id and condition.
    dfa_values{k,1} = png_filename;
    dfa_values{k,2} = a;
    
    close all % Close figures before looping through to the next trial
end

% Exporting DFA values into a table in long format for statistical
% analysis in R.
dfa_results = cell2table(dfa_values, 'VariableNames', {'id','alpha'});
writetable(dfa_results, fullfile(output_directory, 'COP DFA Values.csv'));

% Done!

