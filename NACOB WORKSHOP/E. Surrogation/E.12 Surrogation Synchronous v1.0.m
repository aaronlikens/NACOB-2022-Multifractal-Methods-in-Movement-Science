% E.12 Surrogation Synchronous
% Description
% - This m-file accompanies the synchonous problem set for the same topic.
%   While it does not provide the answers persay it does provide code to
%   help answer them
%% 1 Basic statistics

% This code demonstrates some basic statistics with MATLAB that can be used
% with surrogation. The code can be easily modified and run repetitively to
% get a feel for how this all works.

clear, clc

% This value sets the baseline difference between our original statistic
% and surrogate statistics.
dx=0;

% This value will change the distribution of our surrogate statistics.
A=0.5;

% This is our pretent original statistic. When this array is expanded below
% the original statistic will be the first index.
SE=1.2;

% This step performs the pretend surrogation by creating random numbers.
SE(2:21,1)=SE(1)+dx+A*randn(20,1);

% Performs the rank order.
[SE_s,ind]=sort(SE);

% Tests the null hypothesis.
[~,p]=ttest(SE(2:end),SE(1));

% This plot will specifically make the original statistic a different
% color.
figure
plot(1:21,SE_s,'.b',find(ind==1,1),SE(1),'r.','MarkerSize',15)
xlabel('Rank Order')
ylabel('Descriminating Statistic')
title({'Results of a Pretend Surrogate Test',['p = ' num2str(p)]})
axis tight
xticks(1:21)

%% 2 Thieler surrogation

% Here we will perform Thieler surrogation for one of the algorithms. The
% algorithm used can be changed easily and the plots remade.

clear, clc

% Noise level in the linear stochastic process.
noise=1000;
% Weight in the Linear stochastic process.
weight=1;
% Create a linear stochastic series
y=zeros(2000,1);
for i=2:2000
    y(i)=weight*y(i-1)+noise*randn(1);
end

% Finds the Sample Entropy of the original data.
m=2;
R=0.2;
SE(1)=Ent_Samp20200723(y,m,R);

% We'll use algorithm 0.
algorithm = 0;
% Finds the surrogates and Sample Entropy of those surrogates.
for i=1:20
    y_surr=Surr_Theiler20200723(y,algorithm);
    SE(i+1,1)=Ent_Samp20200723(y_surr,m,R);
end

% Performs the rank order.
[SE_s,ind]=sort(SE);

% Tests the null hypothesis.
[~,p]=ttest(SE(2:end),SE(1));

% Create the figure.

figure

subplot(4,1,1)
plot(y)
title('Original')

subplot(4,1,2)
plot(y_surr)
title('Example Surrogate')

subplot(4,1,[3,4])
plot(1:21,SE_s,'.b',find(ind==1,1),SE(1),'r.','MarkerSize',12)
xlabel('Rank Order')
ylabel('Sample Entropy')
title({['Results of Surrogate Analysis with Algorithm ' num2str(algorithm)],['p = ' num2str(p)]})
axis tight
xticks(1:21)

%% 3 Pseudoperiodic surrogation

% This code will create and display some Pseudoperiodic surrogates.

clear, clc

% This data is an ankle angle during walking.
load('E.12 Surrogation Synchronous v1.0.mat')

% We've hard-coded these in so the AMI and FNN subroutines are not needed.
tau=15;
dim=6;

% Find the optimal rho value.
[rho,out]=Surr_findrho20200630(x,tau,dim);

% Create a figure and the surrogates
figure
subplot(5,1,1),plot(x)
ylabel('O','FontSize',15)

for i=1:4
    % Perform the surrogation
    [xs,xi]=Surr_PseudoPeriodic20200630(x,tau,dim,rho);
    subplot(5,1,i+1),plot(xs)
    ylabel(['PPS' num2str(i)],'FontSize',15)
    axis tight
    
end
