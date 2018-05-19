function slq_algo1
%% writing a good code
% all functions, structs and classes- Camel case
% all variables-underscore
% oldpath=path;
% path(oldpath,'myAutomaticDifferentiation')

%% main
modelParams=setParams();
% load initial trajectory
xInit=[0;0]
load('trajectory.mat', 'trajectory');
trajectory.x=trajectory.x';
trajectory.x=[xInit trajectory.x];
trajectory.u=trajectory.u';

%repeat until max number of iterations or converged (l(t)<l_t)
max_iter=9999;
while max_iter<10000
    %simulate the trajectory
    for time_iter=1:1:modelParams.T/modelParams.dt
        [~, xNext(:,time_iter)]=simplePendDynamics(xInit,trajectory.u(:,time_iter),modelParams);
        xInit=xNext(:,time_iter);
    end
    trajectory.x=[xInit xNext];
    % linearize the dynamics
    [A,B]=linDynamics(modelParams,trajectory);
    % Quadratize cost function along the trajectory
    %--is it optional?----
    %Backwards solve the Ricatti-like equations
    P=zeros(2,2,modelParams.N+1);
    P(:,:,modelParams.N+1)=modelParams.Qf;
    for ricatti_iter=modelParams.N:-1:1
        H=modelParams.Rt+B(:,ricatti_iter)'*P(:,:,ricatti_iter+1)*B(:,ricatti_iter);
        G=B(:,ricatti_iter)'*P(:,:,ricatti_iter+1)*A(:,:,ricatti_iter);
        K(ricatti_iter,:)=-inv(H)*G;
        P(:,:,ricatti_iter)=modelParams.Qt+...
            A(:,:,ricatti_iter)'*P(:,:,ricatti_iter+1)*A(:,:,ricatti_iter)...
            +K(ricatti_iter,:)'*H*K(ricatti_iter,:)+...
            +K(ricatti_iter,:)'*G+ G'*K(ricatti_iter,:);
%         p(:,ricatti_iter)=q;    
    if norm(l)<10e-5
        break
    end
    max_iter=max_iter+1;
end

%% Given
    function modelParams=setParams()
        modelParams.g=1;%9.8
        modelParams.length=1;
        modelParams.dt=0.1;
        modelParams.T=10; %N=T/dt
        modelParams.N=modelParams.T/modelParams.dt;
        modelParams.Qt=diag([10,10]);
        modelParams.Qf=diag([100,100]);
        modelParams.Rt=1;
    end

    %% dynamics of a simple pendulum
    % params x: state vector 2x1
    % params u: input vector 1
    % params modelParams: struct
    %returns xNext: 2x1, xdot=1x2
    function [xdot, xNext]=simplePendDynamics(x,u, modelParams)
        xdot(1)=x(2);
        xdot(2)=-(modelParams.g/modelParams.length)*sin(x(1))+u;
        xNext=x+xdot'*modelParams.dt;
    end

    %% cost function : objective is to minimize this
    % params x: state vector 2x(N+1) -> (0,dt,...,T) where N=T/dt
    % params u: input vector 1x(N)
    % params trajectory: struct (x,u)
    function J = costFunction(x,u, modelParams, trajectory)
        recurr_cost=0;
        for t=1:size(x,2)-1
            error_traj=x(:,t)-trajectory.x(:,t);
            recurr_cost=recurr_cost+ error_traj'*modelParams.Qt*error_traj+...
                u(:,t)'*modelParams.Rt*u(:,t);
        end
        error_goal=x(:,end)-trajectory.x(:,end);
        J= error_goal'*modelParams.Qf*error_goal + recurr_cost;
    end
   
    %% linearization of system dynamics at trajectory points
    %returns A 2x2xN B 2xN
    function [A,B]=linDynamics(modelParams,trajectory)
        delta=0.01;
        % calculate A
        for traj_iter=1:1:(modelParams.T/modelParams.dt)
            [f_x1_plus,~]=simplePendDynamics(trajectory.x(:,traj_iter)+[delta;0],...
                trajectory.u(:,traj_iter),modelParams);
            [f_x1_min,~]= simplePendDynamics(trajectory.x(:,traj_iter)-[delta;0],...
                trajectory.u(:,traj_iter), modelParams);
            [f_x2_plus,~]=simplePendDynamics(trajectory.x(:,traj_iter)+[0;delta],...
                trajectory.u(:,traj_iter),modelParams);
            [f_x2_min,~]= simplePendDynamics(trajectory.x(:,traj_iter)-[0;delta],...
                trajectory.u(:,traj_iter), modelParams);
            [f_u_plus,~]=simplePendDynamics(trajectory.x(:,traj_iter),...
                trajectory.u(:,traj_iter)+delta,modelParams);
            [f_u_min,~]= simplePendDynamics(trajectory.x(:,traj_iter),...
                trajectory.u(:,traj_iter)-delta, modelParams);
            B(:,traj_iter)=(f_u_plus'-f_u_min')/(2*delta);
            A(:,:,traj_iter)=[(f_x1_plus'-f_x1_min')/(2*delta) ...
                (f_x2_plus'-f_x2_min')/(2*delta)];
        end        
    end
end