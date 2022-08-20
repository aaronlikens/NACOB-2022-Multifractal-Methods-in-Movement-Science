
%% 
% DESCRIPTION: 
% This file accompanies the slide deck named "J.01 Multifractal Analysis v1.0 Part 3". This code will make
% references to slides from this presentation. 

% NOTE: This codes uses a function that is not currently in the Nonlinear
% Analysis Core GitHub repository. 

% Section I: Tutorial Example
%             - Contains step by step infomation on how to run MFDFA Analysis.

% Section II: Batch Processing Example
%             - Load all 3 subjects
%             - Run the analysis for all subjects using a for loop
%             - Create and save a summary table

% Required Toolboxes:
% Statistics and Machine Learning Toolbox
% Signal Processing Toolbox
% Image Processing Toolbox


%% SECTION 1 - TUTORIAL EXAMPLE

% -------------------- Step 1:Load Stride Interval ------------------------
%add folder to the path****
% Select location of stored data. In this case, it is in the "DATA" folder.
% We have used 'uigetdir' here as paths for everyone are different.
my_directory = uigetdir(matlabroot, 'select data location');

% Source of data: Raffalt et al., 2021
load S206_selfPaced_StrideIntervals.mat

% Put the time and intervals into their own variables.
time = SI(:,1);
StrideIntervals = SI(:,2);

% ------------ Step 2: Visual inspection before analysis ------------------

% It is always a good idea to visually inspect your data before running any
% analysis and especially when doing pilot collections. This will help to
% make sure you identify any abnormalities and avoid time and money loss.

%   What to look for in this data?
%   At this point of the pre-processing, your data should look like
%   a "noisy" time series (refer to slide 25-26 of the
%  Multifractal presentation for best practices).


figure(1);
plot(StrideIntervals, 'r-')
title('Stride Intervals from Overground Walking Trial', 'FontSize', 16)
xlabel('Time (s)');
ylabel('Stride Interval (s)');
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
disp('Press any key to continue!')

% ----------- Step 3: Set MFDFA parameters for Analysis ------------------

% Pre-define the inout parameters. Details can be found in the MFDFA1 function
scale = [16, 32, 64, 128, 256]; % Note the increase by power of 2

% DFA looks at q = 2, with MFDFA we can look at negative q-orders giving more
% information about the structure of the data.

q2 = -2:0.1:2; % Second statistical moment
q3 = -3:0.1:3; % Third statistical moment
q5 = -5:0.1:5; % Fifth statistical moment

% Polynomial order for detrending
m = 1; 

% Flag for output plot (1 = on, 0 = off)
plotOption = 1; 

% The next line of code runs the actual MFDFA function. An explanation of
% the steps in the analysis can be seen in slides 11-19 of "Mulitfractals - Part 3". 

figure(2);
[Hq, tq, hq, Dq, Fq] = MFDFA1(StrideIntervals, scale,q3 , m, plotOption);

% In this example, we are particularly interested in the spectrum width. 
% Larger width values can be interpreted as greater number of patterns in 
% the time series.
% The width will increase as the q-order increases, but there are limitations.
% In order to perform MFDFA with higher q-order you need a lot of data points.
% For short time series, a q-order that ranges between -3 to 3 should be
% enough.

% Calculate the MFDFA width
width = max(hq)- min(hq);

%% SECTION 2: BATCH PROCESSING EXAMPLE

% DATA AND EXPERIMENT DESCRIPTION
% Participants were instructed to synchronize their right step to various
% metronomes: Pink noise and White noise. Prior to those trials, participants
% performed a self-paced trial for baseline.

clear;
clc;
close all;
    
% Select location of stored data. In this case, it is in the "DATA" folder.
% We have used 'uigetdir' here as paths for everyone are different.
my_directory = uigetdir(matlabroot, 'select data location');
addpath(my_directory); 

% Returns a structure of the path specified in base path. If there are
% other files types (outputs from multiple pieces of equipment) we can
% specify the data we want to load using the '*' and the file type. This
% will tell MATLAB to look for all files with the extension .mat regardless
% of what they are named. 
my_files = dir(fullfile(my_directory, '*.mat')); 

% We also want to save our output so we can go back and view figures and results
% later without having to run the code again. We will create a folder named
% "ANALYSIS OUTPUT" and store everything in here. 
newStr = erase(my_directory, "/DATA")
output_directory = append(newStr, "/ANALYSIS OUTPUT");

if ~exist(output_directory, 'dir')
    mkdir(output_directory)
end

% Create variables for summary table

id = []; % subject ID
cond = []; % condition name
width =[]; % spectrum width
Hurst_exp = []; % Hurst exponent
mass_exp = []; % Mass exponent

% The easiest way to process many files is with a for loop. For this to
% work it is important that you use the same naming convention for all
% files (it is good practice to do this regardless of how you analyze your
% data).

for p = 1:numel(my_files) 
    
    base_filename = my_files(p).name;
    path_file = fullfile(my_files(p).folder, base_filename);
    
    % Read matlab file
    data = load(path_file); 
    
    % We can use the file name to help label our figures and fill in our
    % results table. 
    name = strsplit(base_filename, '_');
    ID = name{1,1};
    Cond = name{1,2};
    
    % Display ID and condition to keep track of the progress of the loop
    disp(append(ID, " -- ", Cond)); 
    
    % Create a time variable: select the first column of the stride intervals
    % variable (SI).
    time = data.SI(:,1);
    
    % Select the second column of the stride intervals
    % variable (SI) to create a 'StrideIntervals' variable.
    StrideIntervals= data.SI(:,2);
    
    % Visual inspection before analysis
    % Plot and Inspect your data:
    %      It's important that you look at your graph before analysis.
    %      Make sure that you don't see any abnormal large peaks or other
    %      abnormalities.
    
    %   What to look for?
    %      At this point of the pre-processing, your data should look
    %      noise-like time series (refer to slide 25-26 of the
    %      power-point for best practices.
    
    figure;
    plot(StrideIntervals, 'r-')
    title('Stride Intervals from Overground Walking Trial', 'FontSize', 16)
    xlabel('Time (s)');
    ylabel('Time Interval (s)');
    pause(2);
    close all; 
    
    % Set MFDFA Parameters for Analysis
    scale = [16, 32, 64, 128, 256];
    
    % DFA looks at q = 2, with MFDFA we can look at negative q-order giving more
    % information about the structure.
    
    q2 = -2:0.1:2; % Second statistical moment
    q3 = -3:0.1:3; % Third statistical moment
    q5 = -5:0.1:5; % Fifth statistical moment
    
    % Polynomial order for detrending
    m = 1;
    
    % Flag for output plot (1 = on, 0 = off)
    plotOption = 1;
    
    % MFDFA starts here
    [Hq,tq,hq,Dq,Fq] = MFDFA1(StrideIntervals, scale,q5 , m, plotOption);
    set(gcf, 'Position', [100, 100, 1000, 800]);
    
    % Export and save our figures into the folder titled "ANALYSIS OUTPUT"
    image = fullfile(output_directory, append(ID, Cond, '.png'));
    saveas(gcf, image);

    pause(2);
    close all;
    
    % Hq: General Hurst exponent
    % tq: Mass exponent
    MFDFA_width= max(hq)- min(hq);
    Hq_avg = mean(Hq);
    tq_avg = mean(tq);
    
    
    % Append all the variables in one table
    id = [id; string(ID),];
    cond = [cond; string(Cond),];
    width = [width; MFDFA_width,];
    Hurst_exp= [Hurst_exp; Hq_avg,];
    mass_exp= [mass_exp; tq_avg,];

end

% SUMMARY TABLE
T = table(id, cond, width, Hurst_exp, mass_exp);
% T = table(ID_cond, width, Hurst_exp, mass_exp);

% Save subject results
writetable(T,append(output_directory, '/MFDFA_ResultsSummary.csv'))

