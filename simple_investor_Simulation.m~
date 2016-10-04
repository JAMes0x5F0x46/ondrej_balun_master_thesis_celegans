% Simulation of Tap Withdrawal (TW) circuit  of C Elegans
% Wicks et al. model
% ID of Neurons: 
    % MEM: 1
    % MRK: 2
    % CBUY: 3
    % CSHT: 4
    % BUY: 5
    % SHT: 6
% Touch sensory neurons: PLM, ALM and AVM. Stimulus is applied only on
% these neurons.  
function [output] = simple_investor_Simulation(stimID, knockedOut)
%Input: 
% KnockedOut: The ID(s) of the neurons which are to be knockedout
% stimID: The ID(s) of the neurons to which stimulus will be applied
%% Declare global variable
global IStim;
global W_gap;
global W_syn;
global Cap;
global E_syn;
global Res;
global V_eq;
global dt;
% 1 - current value, 2 - historical value, 13 - short, 15 - buy
%% Membrane Properties
Cap = [5;5;15;15;15;15;15;15;5;9.1;9.1;14;15;16;14;16; 15]*1e-12;
% Resistance
Res = [30;30;10;10;10;10;10;10;30;16;16;11;10;9.4;11; 9.4; 10]*1e9;
% Reverse Potential of Neuron based on Polarities
E_syn = [-0.048;-0.048;-0.048;-0.048;0;-0.048;0;-0.048;-0.048;-0.048;-0.048;0;-0.048;0;-0.048;-0.048;-0.048];
number_of_neurons = length(Cap);
%% Loading initial state from file
load 'InitialState_simple_investor2.mat'; % Load intial state in V_init
    %% Network Connectivity Matrices of TW Circuit (Adjacency matrix)
load 'Connectivity_matrix_simple_investor2.mat'; % Will Load W_GAP and W_SYN
%% Adjust Network connectivity for Knockout 
    if knockedOut == 0 % no neuron is knockedout
        W_syn = W_SYN;
        W_gap = W_GAP;
    else % delete the neurons from the graph that are to be knocked out
        % Removing corresponding row and column for adjacency matrix for 
        % Synapse
        W_syn = W_SYN(~ismember(1:number_of_neurons,knockedOut),:); % deleting row first
        W_syn = W_syn(:,~ismember(1:number_of_neurons,knockedOut)); % deleting column 

        % Gap-junction
        W_gap = W_GAP(~ismember(1:number_of_neurons,knockedOut),:); % deleting row first
        W_gap = W_gap(:,~ismember(1:number_of_neurons,knockedOut)); % deleting column 

        % Deleting other Parameters 
        Cap = Cap(~ismember(1:number_of_neurons,knockedOut));
        Res = Res(~ismember(1:number_of_neurons,knockedOut));
        E_syn = E_syn(~ismember(1:number_of_neurons,knockedOut));
        V_init = V_init(~ismember(1:number_of_neurons,knockedOut));
    end

%% Computing Equillibrium potential of each neurons
% Parameter
% Leakage voltage
V_leak = -0.035;
% Gap-junction parameter
g_gap = 5e-9;
%Synaptic Parameters
g_syn = 6e-10; 
N = size(W_syn,1); % Number of Neuron in TW circuit after knockout
% Construct Matrix A
A = zeros(N,N);
B = zeros(N,1);
for i=1:N
    for j=1:N
        if i~=j
            A(i,j) = -Res(i)*W_gap(i,j)*g_gap;
        else
            A(i,j) = 1+Res(i)*sum(W_gap(i,:)*g_gap+0.5*W_syn(i,:)*g_syn);
        end
    end
    B(i) = V_leak+Res(i)*sum(W_syn(i,:)*E_syn*g_syn*0.5);
end
V_eq = linsolve(A,B);

%% Simulation of TW circuit 
simdur = 0.4; % 80ms
dt = 0.0001; % 0.1ms

startTime = 0.00; % time when stimulation is started  10ms
% pulse_dur = 0.15; % stimulation duration 40ms

% Constructing Stimulus Input
simStep = simdur/dt;
IStim = zeros(N,simStep+1);

startInd1 = 1; 
endInd1 = startInd1 + 500;

startInd2 = endInd1;
endInd2 = startInd2 + 500;

startInd3 = endInd2;
endInd3 = startInd3 + 500;

startInd4 = endInd3;
endInd4 = startInd4 + 500;

startInd5 = endInd4;
endInd5 = startInd5 + 500;

startInd6 = endInd5;
endInd6 = startInd6 + 500;

startInd7 = endInd6;
endInd7 = startInd7 + 500;

pulse_stim1 = 2e-10; % Stimulation strength
pulse_stim2 = 8e-10; % Stimulation strength
pulse_stim3 = 4e-10; % Stimulation strength
pulse_stim4 = 1e-10; % Stimulation strength
pulse_stim5 = 10e-10; % Stimulation strength
pulse_stim6 = 5e-10; % Stimulation strength
pulse_stim7 = 7e-10; % Stimulation strength


IStim(1,startInd1:endInd1) = pulse_stim1;
IStim(2,startInd1:endInd1) = 0;

IStim(1,startInd2:endInd2) = pulse_stim2;
IStim(2,startInd2:endInd2) = pulse_stim1;

IStim(1,startInd3:endInd3) = pulse_stim3;
IStim(2,startInd3:endInd3) = pulse_stim2;

IStim(1,startInd4:endInd4) = pulse_stim4;
IStim(2,startInd4:endInd4) = pulse_stim3;

IStim(1,startInd5:endInd5) = pulse_stim5;
IStim(2,startInd5:endInd5) = pulse_stim4;

IStim(1,startInd6:endInd6) = pulse_stim6;
IStim(2,startInd6:endInd6) = pulse_stim5;

IStim(1,startInd7:endInd7) = pulse_stim7;
IStim(2,startInd7:endInd7) = pulse_stim6;

%% Call to Simulation method (ODE45 method)
tic;
[store_t,store_V] = ode45('TWModel_dynamics', 0:dt:simdur, V_init);
t_ode = toc;

%% plotting figures

h=figure; hold on;
subplot(1,1,1); hold on;
% plot(store_t*1000, store_V(:,1)*1000,'black','linewidth', 2);
% plot(store_t*1000, store_V(:,2)*1000,'g','linewidth', 2);
plot(store_t*1000, store_V(:,13)*1000,'y','linewidth', 2);
plot(store_t*1000, store_V(:,15)*1000,'c','linewidth', 2);

legend('current', 'historical', 'SHORT', 'BUY');
xlabel('Time (ms)');
ylabel('Potential (mV)');
title('Trajectories of membrane potentials');
% 
% disp(store_V(2001,1));
% disp(store_V(2001,2));
% disp(store_V(2001,3));
% disp(store_V(2001,4));
% disp(store_V(2001,5));
% disp(store_V(2001,6));
% disp(store_V(2001,7));
% disp(store_V(2001,8));
% disp(store_V(2001,9));
% disp(store_V(2001,10));
% disp(store_V(2001,11));
% disp(store_V(2001,12));
% disp(store_V(2001,13));
% disp(store_V(2001,14));
% disp(store_V(2001,15));
% disp(store_V(2001,16));
% disp(store_V(2001,17));
output = [store_t * 1000, store_V(:,1)*1000, store_V(:,2)*1000, store_V(:,13)*1000, store_V(:,15)*1000];

end