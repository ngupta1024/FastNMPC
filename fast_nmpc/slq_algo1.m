function slq_algo1
%% writing a good code
% all functions, structs and classes- Camel case
% all variables-underscore
close all
%% main
modelParams=setParams();

%% handles for dynamics
simple_pend=@(x,u)simplePendDynamics(x,u, modelParams);
aug_pend=@(x,u)augPendDynamics(x,u, modelParams);

%% generate trajectory
%to get the full trajectory- add x_init in each traj.x
if modelParams.gen_traj
    [nom_traj.x, nom_traj.u]=generateTraj(modelParams, simple_pend,[pi;0] );
    save('nominal.mat','nom_traj');
    [des_traj.x,des_traj.u]=generateTraj(modelParams, aug_pend, [3*pi/4;0.2]);
    save('desired.mat','des_traj');
else
    load('nominal.mat','nom_traj');
    load('desired.mat','des_traj');
end
if modelParams.viz
    set(0,'DefaultFigureWindowStyle','docked')
    plot([1:1:modelParams.N],nom_traj.u)
    hold on;
    plot([1:1:modelParams.N],des_traj.u)
    legend('nominal','desired');
    hold off;
    figure;
    plot(nom_traj.x(1,:),nom_traj.x(2,:))
    hold on;
    plot(des_traj.x(1,:),des_traj.x(2,:))
    legend('nominal','desired');
    hold off;
    
end
%% loop
%repeat until max number of iterations or converged (l(t)<l_t)
max_iter=9999;
while max_iter<10000
    %simulate the trajectory
    xInit=modelParams.x_init;
    for time_iter=1:modelParams.T/modelParams.dt
        [~, xNext(:,time_iter)]=simplePendDynamics(xInit,nom_traj.u(:,time_iter),modelParams);
        xInit=xNext(:,time_iter);
    end
    nom_traj.x=[modelParams.x_init xNext];
    des_traj.x=[modelParams.x_init des_traj.x];
    
    % linearize the dynamics
    [A,B]=linDynamics(modelParams,nom_traj);
    
    % Quadratize cost function along the trajectory
    p1=zeros(2,2,modelParams.N+1);
    p1(:,:,modelParams.N+1)=modelParams.Qf;
    p2=zeros(2,modelParams.N+1);
    p2(:,modelParams.N+1)=2*modelParams.Qf*(nom_traj.x(:,end)-des_traj.x(:,end));
    x_diff=nom_traj.x-des_traj.x;
    u_diff=nom_traj.u-des_traj.u;
    q_t=2*bsxfun(@times, diag(modelParams.Qt), x_diff);
    r_t=2*bsxfun(@times, diag(modelParams.Rt), u_diff);
    
    %initialize ricatti variables
%     G=zeros(modelParams.N,2);
%     g=zeros(modelParams.N);
    K=zeros(modelParams.N,2);
    l=zeros(modelParams.N,1);
    
    %Backwards solve the Ricatti-like equations
    for ricatti_iter=modelParams.N:-1:1
        %H-scalar
        H=modelParams.Rt+B(:,ricatti_iter)'*p1(:,:,ricatti_iter+1)*B(:,ricatti_iter);
        %G-1x2
        G=B(:,ricatti_iter)'*p1(:,:,ricatti_iter+1)*A(:,:,ricatti_iter);
        %g-scalar
        g=r_t(ricatti_iter)+B(:,ricatti_iter)'*p2(:,ricatti_iter+1);
        %K-1x2
        K(ricatti_iter,:)=-inv(H)*G;
        %l-scalar
        l(ricatti_iter)=-inv(H)*g;
        %p1-2x2
        p1(:,:,ricatti_iter)=modelParams.Qt+...
            A(:,:,ricatti_iter)'*p1(:,:,ricatti_iter+1)*A(:,:,ricatti_iter)...
            +K(ricatti_iter,:)'*H*K(ricatti_iter,:)+...
            +K(ricatti_iter,:)'*G+ G'*K(ricatti_iter,:);
        %p2-2x1
        p2(:,ricatti_iter)=q_t(:,ricatti_iter)+...
            A(:,:,ricatti_iter)'*p2(:,ricatti_iter+1)...
            + K(ricatti_iter,:)'*H*l(ricatti_iter)+...
            l(ricatti_iter)'*g+...
            G'*l(ricatti_iter);
        p1(:,:,ricatti_iter)
        p2(:,ricatti_iter)
    end
    max_iter=max_iter+1;
end

    %% cost function : objective is to minimize this
    % params x: state vector 2x(N+1) -> (0,dt,...,T) where N=T/dt
    % params u: input vector 1x(N)
    % params trajectory: struct (x,u)
    % function J = costFunction(x,u, modelParams, trajectory)
    %     recurr_cost=0;
    %     for t=1:size(x,2)-1
    %         error_traj=x(:,t)-trajectory.x(:,t);
    %         recurr_cost=recurr_cost+ error_traj'*modelParams.Qt*error_traj+...
    %             u(:,t)'*modelParams.Rt*u(:,t);
    %     end
    %     error_goal=x(:,end)-trajectory.x(:,end);
    %     J= error_goal'*modelParams.Qf*error_goal + recurr_cost;
    % end
    
    
end