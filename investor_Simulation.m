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
function [output] = investor_Simulation(stimID, knockedOut)
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

%% Membrane Properties
% Capacittance
Cap = [5;9.1;14;16;15;16;5;5;5;5]*1e-12;
% Resistance
Res = [30;16;11;9.4;10;11;30;30;30;30]*1e9;
% Reverse Potential of Neuron based on Polarities
E_syn = [-0.048;-0.048;0;0;-0.048;-0.048;-0.048;-0.048;-0.048;-0.048];
number_of_neurons = length(Cap);
%% Loading initial state from file
load 'InitialState_investor.mat'; % Load intial state in V_init
    %% Network Connectivity Matrices of TW Circuit (Adjacency matrix)
load 'Connectivity_matrix_investor.mat'; % Will Load W_GAP and W_SYN
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
simdur = 0.2; % 80ms
dt = 0.0001; % 0.1ms

startTime = 0.01; % time when stimulation is started  10ms
pulse_dur = 0.15; % stimulation duration 40ms

% Constructing Stimulus Input
simStep = simdur/dt;
IStim = zeros(N,simStep+1);

startInd = startTime/dt; 
endInd = (startTime+pulse_dur)/dt;

endInd1 = startInd + 400;
startInd1 = endInd1 + 150;
endInd2 = startInd1 + 400;
startInd2 = endInd2 + 150;

pulse_stim1 = 30e-10; % Stimulation strength
pulse_stim2 = 3e-10; % Stimulation strength
pulse_stim3 = 17e-10; % Stimulation strength
%IStim(stimID,sta2tInd:endInd) = pulse_stim;
%pulse_stim = (rand(1, length(IStim(stimID,startInd:endInd))) * pulse_stim);
IStim(stimID,startInd:endInd1) = pulse_stim1;
IStim(stimID,startInd1:endInd2) = pulse_stim2;
IStim(stimID,startInd2:endInd) = pulse_stim3;

%% Call to Simulation method (ODE45 method)
tic;
[store_t,store_V] = ode45('TWModel_dynamics', 0:dt:simdur, V_init);
t_ode = toc;

%% Getting Circuit Output
% Determine new ID of AVA and AVB (may change for knockout)
%[AVA_ID, AVB_ID] = newNeuronID(knockedOut);
% output = sum_{stim duration}(V_{AVB}-V_{AVA})
% Output > 0 -> Forward/Accelerate movement
% Output < 0 -> Backward movement

output_a = sum((store_V(startInd:endInd,5) - ...
    store_V(startInd:endInd,6)))*(startTime+pulse_dur);

output_buy = store_V(:,5);
output_sht = store_V(:,6);
% disp(output_buy(800));
% disp(output_sht(800));
% disp(store_V(800,2));
% disp(store_V(800,3));
% disp(store_V(800,4));

forward = 0;
backward= 0;
%for i=1:length(output_avb) %entire simulation
for i=100:400 %only stimulation
    if output_buy(i) >= output_sht(i)
        forward = forward + 1;
    else
        backward = backward + 1;
    end
end
output_f = forward; 
output_b= backward; 

%% plotting figures

% Plotting V_{AVA} and V_{AVB}, time and voltage is multiplied by 1000 to 
% convert mS and mV unit respectively
% 
h=figure; hold on;
subplot(1,1,1); hold on;
plot(store_t*1000, store_V(:,5)*1000,'b','linewidth',2); 
plot(store_t*1000, store_V(:,6)*1000,'r','linewidth', 2);
%plot(store_t*1000, store_V(:,2)*1000,'g','linewidth', 2);
%plot(store_t*1000, store_V(:,3)*1000,'y','linewidth', 2);
%plot(store_t*1000, store_V(:,4)*1000,'c','linewidth', 2);
%plot(store_t*1000, store_V(:,1)*1000,'black','linewidth', 2);
%plot(store_t*1000, store_V(:,7)*1000,'black','linewidth', 2);
%plot(store_t*1000, store_V(:,10)*1000,'black','linewidth', 2);
legend('BUY','SHORT', 'MRK', 'CBUY', 'CSHT');
xlabel('Time (ms)');
ylabel('Potential (mV)');
title('Trajectories of membrane potentials');
% 
% % plotting stimulus
% subplot(2,2,3);
% plot(store_t*1000, IStim(stimID(1),:),'linewidth',2);
% xlabel('Time (ms)');
% ylabel('Current (Amp)');
% title('Stimulus Current');
% fname = ['Output_',num2str(stimID),'_',num2str(knockedOut),'.eps'];


output = [output_f, output_b];
%output = [store_t, store_V(:,AVA_ID) * 1000, store_V(:,AVB_ID) * 1000];

end