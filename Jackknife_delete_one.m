%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is to replicate paper "The Jackknife-a review"
% Author: Guangyao Zhou & Yutao Sun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replication results are stored in some Excel files. 
function result = Jackknife_delete_one
    for j = 1:19
       file_name = "Result_Jackknife_" + string(j);
       result = repetition(j);
       % "j" determines the value of mu2. In this study, mu2 have 19 values and
       % denotes the 
       xlswrite(file_name,result,"Sheet1","A1");
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function repl = repetition(j)
% In our study, T is the simulation times; "table" consists of 2SLS
% estimator, Jackknife-delete-one estimator, and Coverage of boostrapping
% interval;
T = 1000;
table = cell(T,3);
for i = 1:T % The parfor code here is parfor i = 1:T
    result_temp = main(j);
    table(i,:) = result_temp(:);
end
repl = table;
end

function stats = main(j)
% Main function.Collect the final results.
Data = initialization(j);
Result_2SLS = twostage(Data);
Result_1 = JackknifeDeleteOne(Data);
temp1 = Result_1;
resampling_jack_1 = [];
for i = 1:999
Data_resampling = bootstrapping(Data); % Check 4
Result_Del_2SLS = twostage(Data_resampling);
Result_Del_One_Resampling = JackknifeDeleteOne(Data_resampling);
resampling_jack_1 = [resampling_jack_1;Result_Del_One_Resampling];
end
coverage_one = Coverage(resampling_jack_1);

stats = [{Result_2SLS},{Result_1},{coverage_one}];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2SLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = twostage(Data)
    % Get 2SLS estimator;
    y = Data(:,1);
    x = Data(:,2);
    z = Data(:,3:5);
    pi = z\x; Zpi = z*pi;
    result = Zpi\y;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jackknife_Delelte_One
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = JackknifeDeleteOne(Data)
    sum_jack_one_1 = 0;
    Beta = twostage(Data);
    for i = 1:500
            % Delete ith elements and perform 2SLS model.
            Data_tmp = Data;
            Data_tmp(i,:)=[];
            Jack_one_Beta = twostage(Data_tmp);
            sum_jack_one_1 = sum_jack_one_1 + Jack_one_Beta;
    end
    result = Beta*500 - 499*(sum_jack_one_1/500.0);
    % This result contains two components: The first one is to show the
    % final Jackknife-delete-one estimator; The second one is a template
    % value that is helpful with the calculation of Jackknife-delete-two
    % estimators and Jackknife-delete-three estimators.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrapping 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = bootstrapping(Data)
%replacement
size_Data = size(Data);
n=size_Data(1);
idx= ceil(n*rand(n,1)); %generate n random index between 1 and n
result = Data(idx,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coverage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = Coverage(Data)
% Find the coverage of the given interval; Zero in our case is the true value.
    BS_interval = sort(Data,'ascend');
    if (BS_interval(25,1) <= 0 && BS_interval(975,1) >=0)
        result = 1;
    else
        result = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = initialization(j)
% Parameter : mu2 = [0, 0.25, 1 : 15, 30, 1000]; 
% mu2 indicates the strength of the instruments (Stock, Wright, & Yogo,
% 2002).
% beta = 0;
% Error term~Normal(0, 1);
templ=(1:15)';
allmu2 = [0; 0.25; templ; 30; 1000];
mu2=allmu2(j);
R = mvnrnd([0;0], [1, 0.99; 0.99, 1], 500);
Z = normrnd(0,1,500,3);
PI = ones(3, 1); 
a = sqrt((mu2)/(PI'*Z'*Z*PI));
PI=PI*a;
u=R(:, 1);
v=R(:, 2);
x=Z*PI+v;
beta = 0;
y=beta + u;
result = [y,x,Z];
end
