function slqMPC_algo2
%% writing a good code
% all functions, structs and classes- Camel case
% all variables-underscore
close all
modelParams=setParams();

%% initialization
curr_state= modelParams.x_init;
curr_input=0;
actual_traj.x=zeros(2,modelParams.N);
actual_traj.u=zeros(1,modelParams.N);
N_initial=modelParams.N;

while modelParams.T>1e-5
    %% Adjust Setup
    nom_traj.x=zeros(2,modelParams.N);
    nom_traj.u=zeros(1,modelParams.N);
    % if time lag is there, predict state by forward simulating dynamics with
    % previous controller 

    %% initialize SLQ with controller

    %linearize dynamics about curr_state and curr_input

    curr_traj.x=curr_state;
    curr_traj.u=curr_input;
    [A_0,B_0]=linDynamics(modelParams,curr_traj);

    %compute LQR at linearized state
    [K_0,~]=lqr(A_0, B_0, modelParams.Q_lqr, modelParams.Rt);

    %% solve using SLQ
    
    %forward simulation
    x_desired=modelParams.goal;
    nom_traj.x(:,1)=curr_state;
    for fwd_iter=1:modelParams.N-1
        nom_traj.u(fwd_iter)=-K_0*(nom_traj.x(:,fwd_iter)-x_desired);
        [~,nom_traj.x(:,fwd_iter+1)]=simplePendDynamics(nom_traj.x(:,fwd_iter), nom_traj.u(fwd_iter), modelParams);
    end
    nom_traj.u(modelParams.N)=0;
    
    %run SLQ
    sim_traj=slq_algo1(nom_traj,modelParams);

    % plan for 10 seconds in SLQ and use mpc_steps only
    actual_traj.x(:,N_initial-modelParams.N+1:N_initial-modelParams.N+modelParams.mpc_steps)=sim_traj.x(:,1:modelParams.mpc_steps);
    
    %plot actual trajectory on the desired trajectory plot
    figure(2);
    hold on;
    plot(actual_traj.x(1,1:N_initial-modelParams.N+modelParams.mpc_steps),...
        actual_traj.x(2,1:N_initial-modelParams.N+modelParams.mpc_steps),...
        'c','LineWidth',2,'DisplayName','Actual');
    legend('show');
    drawnow
    
    %% update curr state, time horizon and stuff
    curr_state=sim_traj.x(:,modelParams.mpc_steps+1);
    curr_input=sim_traj.u(:,modelParams.mpc_steps+1);
    modelParams.T=modelParams.T-modelParams.dt*modelParams.mpc_steps;
    modelParams.N=floor(modelParams.T/modelParams.dt)+1;
    modelParams.x_init=curr_state;
end
save('actual_traj.mat','actual_traj');
end