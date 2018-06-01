function slq_algo1
%% writing a good code
% all functions, structs and classes- Camel case
% all variables-underscore
close all

%% main
modelParams=setParams();
last_l=0;
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
    figure(1);
    set(0,'DefaultFigureWindowStyle','docked')
    plot([1:1:modelParams.N],nom_traj.u)
    hold on;
    plot([1:1:modelParams.N],des_traj.u)
    legend('nominal','desired');
    title('control versus time - given')
    hold off;
    figure(2);
    plot(nom_traj.x(1,:),nom_traj.x(2,:))
    hold on;
    plot(des_traj.x(1,:),des_traj.x(2,:))
    legend('nominal','desired');
    title('state space - given')
    hold off;
    drawnow
end
%% loop
%repeat until max number of iterations or converged (l(t)<l_t)
max_iter=1;
while max_iter<10000
    %simulate the trajectory
    nom_traj.x(:,1)=modelParams.x_init;
    for time_iter=1:modelParams.N-1
        [~, nom_traj.x(:,time_iter+1)]=simplePendDynamics(nom_traj.x(:,time_iter),...
            nom_traj.u(:,time_iter),modelParams);
    end
    
    % linearize the dynamics
    [A,B]=linDynamics(modelParams,nom_traj);
    
    %compute cost function
    J_nom=computeActualCost(nom_traj,des_traj,modelParams);
    
    % Quadratize cost function along the trajectory
    p1=zeros(2,2,modelParams.N);
    p1(:,:,modelParams.N)=2*modelParams.Qf;
    p2=zeros(2,modelParams.N);
    p2(:,modelParams.N)=2*modelParams.Qf*(nom_traj.x(:,end)-des_traj.x(:,end));
    x_diff=nom_traj.x-des_traj.x;
    u_diff=nom_traj.u-des_traj.u;
    q_t=2*bsxfun(@times, diag(modelParams.Qt), x_diff);
    r_t=2*bsxfun(@times, diag(modelParams.Rt), u_diff);
    
    %initialize ricatti variables
    K=zeros(modelParams.N-1,2);
    l=zeros(modelParams.N-1,1);
    
    %Backwards solve the Ricatti-like equations
    for ricatti_iter=modelParams.N-1:-1:1
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
            + 2*K(ricatti_iter,:)'*H*l(ricatti_iter)+...
            K(ricatti_iter,:)'*g+...
            2*G'*l(ricatti_iter);
    end
    
    %%Line Search
    alpha=1;
    for ls_iter=1:modelParams.ls_steps
        act_traj.x(:,1)=modelParams.x_init;
        for sim_iter=1:modelParams.N-1
            act_traj.u(sim_iter)= nom_traj.u(sim_iter)+alpha*l(sim_iter)+...
                K(sim_iter,:)*(act_traj.x(:,sim_iter)-nom_traj.x(:,sim_iter));
            if abs(act_traj.u(sim_iter))>modelParams.u_lim
                act_traj.u(sim_iter)=sign(act_traj.u(sim_iter))*modelParams.u_lim;
            end
            [~,act_traj.x(:,sim_iter+1)]=simplePendDynamics(act_traj.x(:,sim_iter),...
                act_traj.u(sim_iter),modelParams);
        end
        act_traj.u(modelParams.N)=0;
        J_actual=computeActualCost(act_traj,des_traj,modelParams);
        if J_actual<J_nom
            break
        end
        alpha=alpha/modelParams.alpha_d;
    end
    if modelParams.viz
        figure(3);
        hold on;
        set(0,'DefaultFigureWindowStyle','docked')
        plot([1:1:modelParams.N],act_traj.u,'DisplayName','actual_'+string(max_iter))
%         hold on;
%         plot([1:1:modelParams.N],nom_traj.u)
        legend('show');
        title('Control Input versus time - actual')
        hold off;
        figure(4);
%         plot(nom_traj.x(1,:),nom_traj.x(2,:))
        hold on;
        plot(act_traj.x(1,:),act_traj.x(2,:),'DisplayName','actual_'+string(max_iter))
        legend('show');
        title('state space - actual')
        hold off;
        drawnow
    end
    nom_traj.u=act_traj.u;
    if abs(last_l-norm(l))<1e-5
        break
    end
    last_l=norm(l)
    max_iter=max_iter+1;
end    
    
end