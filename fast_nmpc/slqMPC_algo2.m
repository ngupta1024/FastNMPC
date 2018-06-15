function slqMPC_algo2
%% writing a good code
% all functions, structs and classes- Camel case
% all variables-underscore
close all
modelParams=setParams();

%% initialization
actTraj.x=modelParams.x_init;
actTraj.u=0;
u_ff=zeros(1,modelParams.N-1);
u_fb=zeros(modelParams.N-1,2);
mpc_iter=0;
% N_initial=modelParams.N;
% time_elapsed=0;

%% loop
while norm(actTraj.x(:,end)-modelParams.goal)>1e-2
    %% Adjust Setup    
    mpc_start=tic;
    % acquire current state 
    curr_state=actTraj.x(:,end);
    curr_input=actTraj.u(:,end);
    modelParams.x_init=curr_state;
    
    % if time lag is there, predict state by forward 
    % simulating dynamics with previous controller 
    
    % update time horizon based on distance
%     modelParams.T= ceil(round(abs(wrapToPi(modelParams.goal(1)-curr_state(1))))*(modelParams.T/abs(wrapToPi(modelParams.x_init(1)-modelParams.goal(1)))));
%     modelParams.N= floor(modelParams.T/modelParams.dt)+1;
    
    %% initialize SLQ with controller

    % linearize dynamics about curr_state and curr_input
    curr_traj.x=curr_state;
    curr_traj.u=curr_input;
    [A_0,B_0]=linDynamics(modelParams,curr_traj);

    % compute LQR at linearized state
    [K_0,~]=lqr(A_0, B_0, modelParams.Q_lqr, modelParams.Rt);

    % initialize SLQ with LQR around current state 
    x_desired=modelParams.goal;
    nom_traj.x(:,1)=curr_state;
    nom_traj.u=zeros(1,modelParams.N);
    for fwd_iter=1:modelParams.N-1
        nom_traj.u(fwd_iter)=-K_0*(nom_traj.x(:,fwd_iter)-x_desired);
        [~,nom_traj.x(:,fwd_iter+1)]=simplePendDynamics(nom_traj.x(:,fwd_iter), nom_traj.u(fwd_iter), modelParams);
    end
    
    %% solve using SLQ

    % run SLQ
    modelParams.viz=0;
    modelParams.printf=0;
    slq_time=tic;
    [nomTraj, u_ff, u_fb]=slq_algo1(nom_traj,modelParams, u_ff, u_fb);
    slq_done=toc(slq_time);
    mpc_done=toc(mpc_start); %equivalent to t_lag
    t_lag=mpc_done;
    time_elapsed=t_lag*2;
    fprintf("the time taken by this MPC iteration = %f for time horizon = %d \n",mpc_done, modelParams.T);
    
    % send control to closed loop controller
    actTrajPart=realPendDynamics(ceil(time_elapsed), modelParams.x_init, nomTraj, u_ff, u_fb, modelParams);
    actTraj.x=[actTraj.x actTrajPart.x];
    actTraj.u=[actTraj.u actTrajPart.u];
    mpc_iter=mpc_iter+1;
end
mpc_iter
fig_pend=figure('Name','slq');
pend_animation(actTraj.x(1,:),fig_pend);
figure('Name','states and inputs');
subplot(2,3,1);
plot([1:length(actTraj.x(1,:))],actTraj.x(1,:),'b','LineWidth',2);
title('position versus time');
subplot(2,3,2);
plot([1:length(actTraj.x(2,:))],actTraj.x(2,:),'b','LineWidth',2);
title('velocity versus time');
subplot(2,3,3);
plot([1:length(actTraj.u)],actTraj.u,'r','LineWidth',2);
title('control input versus time');
subplot(2,3,[4 5 6]);
plot(actTraj.x(1,:),actTraj.x(2,:),'b','LineWidth',2);

% save('actual_traj.mat','actual_traj');
end