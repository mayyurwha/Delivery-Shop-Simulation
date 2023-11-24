clear all; close all;
 
% number of simulations
nsim = 5;
 
% total number of items to be produced before Christmas market  
total_items = 10000;
 
% items in one batch 
batch_size = 100; 
 
% number of batches
n_batches = total_items / batch_size;
 
% 1st line service time (in hours) ~ U([48, 72]) 
a =  48; b = 72; 
a1 =  60; b1 = 80; 
 
% 2nd line service time (in hours) ~ U([72, 120]) 
c = 72; d = 120; 
c1 = 80; d1 = 130;
 
% initial service times for the two service points
service_time = [randi([a, b]) randi([c, d])]; 
service_time1 = [randi([a1, b1]) randi([c1, d1])];
 
% initialize the interarrival time 
inter_arrival = 0; 
 
% time step = 1 hour
timestep = 1; 
 
% If the job is always free to be moved forward
% the corresponding job time is 0
 
% [input, queue, point1, point2, output] for both contractor A&B
job_time = [inter_arrival 0 service_time 0];  %contractor A
job_time1 = [inter_arrival 0 service_time1 0]; %contractor B
 
% number of stages in the system
ns = length(job_time); 
 
% to store number of days needed to complete the order
days_needed = zeros(nsim, 1); %for contractor A
days_needed1 = zeros(nsim, 1); %for contractor B
 
% go through number of simulations
for pp = 1:nsim
    % testing for different batches for contractor A and B
  for jj = 1:n_batches % jj is the number of batches alloted to A
      B_batch = n_batches - jj; %B_batch is the batches alloted to B after deducting jj from total batches
    % initial system state vector for all the stages, 0 idle, 1 occupied
    s = zeros(1,ns);  
    s1 = zeros(1,ns);
    
    % a job waiting in the input
    s (1) = 1; 
    s1 (1) = 1; 
    
    % time the job has been done so far, 
    % equal to 0 for all the stages at the beginning 
    age = zeros(1,ns);
    age1 = zeros(1,ns);
    
    % For monitoring the system
    
    % cumulative output to memory 
    Output = []; % for contractor A
    Output1 = []; % for contractor B
    OutputT = [];  % for both contractors
    
    % State vectors, Ages of the jobs
    S = []; A=[]; % for contractor A
    S1 = []; A1=[]; % for contractor B
    
    % save states of the initial state
    S(1,:) = s; % for contractor A
    S1(1,:) = s1; % for contractor B
    
    % save times done for jobs at the beginning 
    A(1,:) = age;    % for contractor A 
    A1(1,:) = age1;  % for contractor B
    
    % index vector for the stages from which a job can be moved forward
    ii = 1:ns-1;    % for contractor A
    ii1 = 1:ns-1;   % for contractor B
    
    %initialization
    time = 0;
    out_batchesA = 0;
    out_batchesB = 0;
    out_batches = 0;
 
    % while we have some batches to produce
    while out_batches < n_batches 
        
       % states occupied 
       occupied       = s(ii)>0;  % for contractor A
       occupied1       = s1(ii1)>0; % for contractor B
       
       % increase age of ongoing jobs
       age(occupied)  = age(occupied) + timestep;   % for contractor A
       age1(occupied1)  = age1(occupied1) + timestep; % for contractor B
       
       % find which jobs are done
       job_done       = age(ii) >= job_time(ii);   % for contractor A
       job_done1       = age1(ii1) >= job_time1(ii1);  % for contractor B
       
       % find where the next stage is free 
       next_vacant    = s(ii+1)==0;  % for contractor A
       next_vacant1    = s1(ii1+1)==0;  % for contractor B
       
       % queue stage always free for a new job to come!
       next_vacant(1) = 1;    % for contractor A
       next_vacant1(1) = 1;   % for contractor B
       
       % find the indexes of the system for which all these are true               
       move = find(occupied & job_done & next_vacant); % for contractor A
       move1 = find(occupied1 & job_done1 & next_vacant1);  % for contractor B
       
       % do the move in the queue
       
       if length(move)>0 && out_batchesA < jj % or length
           s(move+1)   = s(move+1)+1;  %to this place
           s(move)     = s(move) - 1;  %...from here
           age(move+1) = 0;     % initialize the age counter to zero
           age(move)   = 0;       

           if ismember(3,move+1)  % check for moving to the 1st service line
               job_time(3) = randi([a, b]); % new jobs get new service times
           end
 
           if ismember(4,move+1)  % check for moving to the 2nd service line
               job_time(4) = randi([c, d]);  % new jobs get new service times
           end
 
       end
       
       if length(move1)>0 && out_batchesB < B_batch % or length
           s1(move1+1)   = s1(move1+1)+1;  %to this place
           s1(move1)     = s1(move1) - 1;  %...from here
           age1(move1+1) = 0;     % initialize the age counter to zero
           age1(move1)   = 0;       
 
           if ismember(3,move1+1)  % check for moving to the 1st service line
               job_time1(3) = randi([a1, b1]); % new jobs get new service times
           end
 
           if ismember(4,move1+1)  % check for moving to the 2nd service line
               job_time1(4) = randi([c1, d1]);  % new jobs get new service times
           end
 
       end
      
       % keep input occupied 
       s(1)  = 1;  % for contractor A
       s1(1)  = 1; % for contractor B
 
       % put in memory if a job comes out
       Output= [Output; s(ns)]; % for contractor A
       Output1= [Output1; s1(ns)]; % for contractor B
       OutputT= [OutputT; s(ns)+s1(ns)]; % for both contractors
       
       % keep exit free
       s(ns) = 0;  % for contractor A
       s1(ns) = 0; % for contractor B
       
       %storing of the age and input
       S = [S;s];  % for contractor A
       A = [A;age]; % for contractor A
       S1 = [S1;s1]; % for contractor B
       A1 = [A1;age1]; % for contractor B
 
       % increment time
       time = time + timestep; 
       % calculate how many batches have been done so far
       out_batchesA = sum(Output); % for contractor A
       out_batchesB = sum(Output1); % for contractor B
       out_batches = sum(OutputT);  % for both contractors
    end
 
    % calculate the number of days needed to produce the order
    days_needed(jj, 1) = time / 24; % for contractor A
    days_needed1(jj, 1) = time / 24; % for contractor B
    days_neededT = max(days_needed, days_needed1); % maximum of both contractors, since they are working simultaneously
    
  end
    
    % calculate the average number of days needed to make the order before
    % Christmas market
    
    %computing the job distribution between the contractors based on
    %minimum days required
    min_time_job=find(days_neededT == min(days_neededT,[],'all'));  %job of A for minimum time (batches linked to time jj)
    contractorA_job = min_time_job; % number of batches alloted to A at minimum time
    contractorB_job = n_batches - min_time_job; % % number of batches alloted to B at minimum time
    job_ratio = contractorA_job / contractorB_job; % batch ratio of A to B
    
    days_neededK = days_neededT(min_time_job); % Allocating the first minimum days required for the production since they are the same
    
    days_neededC(pp, 1) = days_neededK(1);   % Accumulated minimum days required for nsim
end

% how many days before Christmas market the contractor should start 
% production to get the order done on time with 95% confidence
days_95 = quantile(days_neededC, 0.95);

days_95_format = ['number of days for 95% confidence is ', num2str(days_95)];

 disp(days_95_format)
 disp('number of jobs for contractor A and B are ')
 disp(contractorA_job * 100)
 disp(contractorB_job * 100)


