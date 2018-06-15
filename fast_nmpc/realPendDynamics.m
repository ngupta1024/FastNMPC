function [actTraj]=realPendDynamics(time, initial_state, nomTraj, u_ff, u_fb, modelParams)
% changing dynamics for real pendulum
modelParams.m=modelParams.m+0.2;
modelParams.g=9.81;
modelParams.length=modelParams.length+0.1;
modelParams.c=modelParams.c;

% simulating dynamics with a closed loop controller
runtime=time/modelParams.dt;
actTraj.u=zeros(1,runtime+1);
actTraj.x=zeros(2,runtime+1);
actTraj.x(:,1)=initial_state;

for fwd_iter=1:runtime
    actTraj.u(fwd_iter)=nomTraj.u(fwd_iter)+u_ff(fwd_iter)+u_fb(fwd_iter,:)*(actTraj.x(:,fwd_iter)-nomTraj.x(:,fwd_iter));
    [~,actTraj.x(:,fwd_iter+1)]=simplePendDynamics(actTraj.x(:,fwd_iter),actTraj.u(fwd_iter),modelParams);
    if abs(actTraj.u(fwd_iter))>modelParams.u_lim
        actTraj.u(fwd_iter)=sign(actTraj.u(fwd_iter))*modelParams.u_lim;
    end
end
actTraj.u(runtime+1)=actTraj.u(runtime);

% fig_pend=figure('Name','slq');
% pend_animation(act_traj.x(1,:),fig_pend);
end