clc;
clear all;
close all;
global b_damp x0 xf N dt nState mControl
global q0 q1 Q P r R q0f q1f Qf S2 S1 S0

global x_nom x_d u_nom u_d t0 uff_col K_col alpha u_temp
nState = 2; % size of state space
mControl = 1; % size of control space

b_damp = 0.3;
N = 100 ;%number of grid points
dt = .1 ;%time step for Euler Integration
T = N*dt ;%total time

% Nominal trajectory
x0 = [0;0];
xf_nom = [pi; 0];
% Desired trajectory
xf_des = [3*pi/4; 0.2];

xf = xf_des;

A = [] ;%empty because no linear equations
b = [] ;%empty because no linear equations
Aeq = [] ;%empty because no linear equations
beq = [] ;%empty because no linear equations
umax = 1 ;%upper control bound (optimization variable) 

options =optimoptions(@fmincon,'TolFun',0.00000001,'MaxIter',10000,...
    'MaxFunEvals',100000,'Display','iter','DiffMinChange',0.001,'Algorithm','sqp');

UB = [inf*ones(1,nState*N),umax*ones(mControl,N)];
LB = -UB;

params0 = [zeros(1,nState*N), zeros(mControl,N)];

t0 = [linspace(0,10,N)]';

% Finding the nominal and desired trajectories
for m = 1:2

    COSTFUN = @(params) cost_to_minimize(params);
    if m==1
        CONSTRAINTFUN = @(params) nonlinconstraints_aug(params,@Aug_Dyn);
    else 
        CONSTRAINTFUN = @(params) nonlinconstraints_aug(params,@Simp_Dyn);
    end
    params = fmincon(COSTFUN,params0,A,b,Aeq,beq,LB,UB,CONSTRAINTFUN,options);

    u = params(2*N+1:3*N);
    x = [params(1:N);params(N+1:2*N)];
    
    U{m} = u;
    X{m} = x;
        
    xf = [pi;0] ;%reset final conditions PART 2b
end

% load X.mat;
% load U.mat;
% 1 is desired, 2 is nominal
%% Plotting desired and nominal trajectories
figure;
for m=1:2
    q = X{m}';
    u = U{m};
    
    subplot(321), hold on %PLOT control input u
       if m==1, plot(t0,u,'LineWidth',2);
       else, plot(t0,u,'LineWidth',2); end

       ylabel({'u(t) '},'fontsize',16)
       title('Control and State Variable Trajectories  ','fontsize',16)
       set(gca,'XTicklabel',[])
    subplot(323), hold on; %PLOT omega response
       if m==1, plot(t0,q(:,2),'LineWidth',2); 
       else, plot(t0,q(:,2),'LineWidth',2); end

       ylabel({'dq(t) '},'fontsize',16)
       set(gca,'XTicklabel',[])
    subplot(325), hold on; %PLOT theta response
       if m==1, plot(t0,q(:,1),'LineWidth',2);
       else, plot(t0,q(:,1),'LineWidth',2); end

       ylabel({'q(t) '},'fontsize',16)
       xlabel('time (sec) ','fontsize',16)
    subplot(3,2,[2,4,6]), hold on; %Phase Plot
       if m==1, plot(q(:,1),q(:,2),'LineWidth',2);
       else, plot(q(:,1),q(:,2),'LineWidth',2); end

       title('Phase Space','fontsize',16)
       xlabel('q(t) ','fontsize',16)
       ylabel('dq(t) ','fontsize',16)
end

%% SLQ Algorithm Buchli
% Global parameters for cost function


% Delaring the global parameters
q0 = 1;
q1 = zeros(nState,1);%ones(nState,1);
Q = eye(nState);
q0f = q0;
q1f = q1;
Qf = 1*Q;

P = zeros(mControl, nState);%ones(mControl, nState);
r = zeros(mControl,1);%ones(mControl,1);
R = 1*eye(mControl);

Alin = zeros(nState,nState,N);
Blin = zeros(nState,mControl,N);
S2=zeros(nState,nState,N);
S1=zeros(nState,N);
S0=zeros(1,N);
outer_iter = 1;
max_iterations_outer = 20;

u_nom = U{2};
x_nom = X{2};

x_d = X{1};
u_d = U{1};

u_temp = [];
figure(2);
   
       subplot(321), hold on %PLOT control input u
       plot(t0,u_nom,'LineWidth',2); hold on;
       plot(t0,u_d,'LineWidth',2);

       ylabel({'u(t) '},'fontsize',16)
       title('Control and State Variable Trajectories  ','fontsize',16)
       set(gca,'XTicklabel',[])
       subplot(323), hold on; %PLOT omega response
       plot(t0,x_nom(2,:),'LineWidth',2); hold on; 
       plot(t0,x_d(2,:),'LineWidth',2);

       ylabel({'dq(t) '},'fontsize',16)
       set(gca,'XTicklabel',[])
    subplot(325), hold on; %PLOT theta response
       plot(t0,x_nom(1,:),'LineWidth',2); hold on;
       plot(t0,x_d(1,:),'LineWidth',2);

       ylabel({'q(t) '},'fontsize',16)
       xlabel('time (sec) ','fontsize',16)
    subplot(3,2,[2,4,6]), hold on; %Phase Plot
       plot(x_nom(1,:),x_nom(2,:),'LineWidth',2); hold on;
       plot(x_d(1,:),x_d(2,:),'LineWidth',2);

       title('Phase Space','fontsize',16)
       xlabel('q(t) ','fontsize',16)
       ylabel('dq(t) ','fontsize',16)
   
while(outer_iter <= max_iterations_outer) % or l(t) <lt

    t0 = linspace(0,10,100)';

%   Linearize A and B about x_nom and u_nom
    
    [Alin, Blin] = linearize();

    % Backward ricatti for S
    [uff_bar, K] = backward_ricatti_forS(Alin, Blin);
    
    uff_col = reshape(uff_bar(:,1),[mControl,1]);
    K_col = reshape(K(:,:,1),[nState*mControl,1]);
    for i = 2:N
       uff_col = [uff_col reshape(uff_bar(:,i),[mControl,1])];
       K_col = [K_col reshape(K(:,:,i),[nState*mControl,1])];
    end
    
%% Forward integrate for u and x - without line search
%     Using linearlized A and B
%     alpha = 1;
%     [x_new, u_new] = forward_integral(uff_bar, K, Alin, Blin, alpha);

%%  Using xdot = f(x)
    alpha = 1;
    [t,x_new] = ode45(@Simp_Dyn_new,t0,x0);
    x_new = x_new';
    
   for i = 1:length(t)
        u_new(:,i) = u_nom(:,i) + alpha*uff_bar(:,i) + K(:,:,i)*(x_new(:,i) - x_nom(:,i));
    end

  %%
    x_nom = x_new;
    u_nom = u_new;
% u_nom = u_temp;
    outer_iter = outer_iter +1;

%     plot(x_nom(1,:),x_nom(2,:),'LineWidth',2); %PLOT x feed forward
%        ylabel({'x(t) '},'fontsize',16)
%        set(gca,'XTicklabel',[])
%     hold on;

       subplot(321), hold on %PLOT control input u
       plot(t0,u_nom,'LineWidth',2); hold on;
       plot(t0,u_d,'LineWidth',2);

       ylabel({'u(t) '},'fontsize',16)
       title('Control and State Variable Trajectories  ','fontsize',16)
       set(gca,'XTicklabel',[])
    subplot(323), hold on; %PLOT omega response
       plot(t0,x_nom(2,:),'LineWidth',2); hold on; 
       plot(t0,x_d(2,:),'LineWidth',2);

       ylabel({'dq(t) '},'fontsize',16)
       set(gca,'XTicklabel',[])
    subplot(325), hold on; %PLOT theta response
       plot(t0,x_nom(1,:),'LineWidth',2); hold on;
       plot(t0,x_d(1,:),'LineWidth',2);

       ylabel({'q(t) '},'fontsize',16)
       xlabel('time (sec) ','fontsize',16)
    subplot(3,2,[2,4,6]), hold on; %Phase Plot
       plot(x_nom(1,:),x_nom(2,:),'LineWidth',2); hold on;
       plot(x_d(1,:),x_d(2,:),'LineWidth',2);

       title('Phase Space','fontsize',16)
       xlabel('q(t) ','fontsize',16)
    pause;
%     pause(0.1);
end
%% NONLINEAR CONSTRAINTS TO MINIMIZE (GENERIC PENDULUM)
function [c,ceq] = nonlinconstraints_aug(params,Dynamics) 
    global N x0 xf dt
    x = [params(1:N);params(N+1:2*N)];
    u = params(2*N+1:3*N);
    
    ceq_0 = x(:,1)-x0;
    ceq_f = x(:,end)-xf;
    ceq = ceq_0;
    x_dot_k = Dynamics(x(:,1),u(1));
    
    for k = 1:N-1        
        x_k = (x(:,k));
        x_k_p1 = x(:,k+1);
        x_dot_k_p1 = Dynamics(x(:,k+1),u(k+1));

        h = dt;
        x_ck = 1/2*(x_k + x_k_p1) + h/8*(x_dot_k - x_dot_k_p1);
        u_ck = (u(k)+u(k+1))/2;
        x_dot_ck = Dynamics(x_ck,u_ck);
      
        defects = (x_k - x_k_p1) + dt/6* (x_dot_k + 4*x_dot_ck + x_dot_k_p1);
        ceq = [ceq;defects];
        x_dot_k = x_dot_k_p1;
    end
    c = [] ;%inequality <= constraint
    ceq = [ceq;ceq_f] ;%equality = constraint
end

%% COST FUNCTION
function J = cost_to_minimize(params)
    global N xf dt %Q R 
    x = [params(1:N);
           params(N+1:2*N)];
    u = params(2*N+1:3*N);

%     Qf = 1*Q;
    J = 0; %(x(:,end)-xf)'*Qf*(x(:,end)-xf);
    Q = eye(2); R = 1;
    for iter = 1:N-1
        J = J+( (x(:,iter))'*Q*(x(:,iter))+u(iter)'*R*u(iter) );%*dt;  
    end
end 

%% SIMPLE Damped PENDULUM DYNAMICS
function dxdt = Simp_Dyn(x,u)
    global b_damp
    dxdt = [x(2);
                u-b_damp*x(2)-sin(x(1))];
end


function dxdt = Aug_Dyn(x,u)
    dxdt = [x(2);
                u-(x(1)^2-1)*x(2)-sin(x(1))];
end

%% SIMPLE Damped PENDULUM DYNAMICS (for ODE)
function dxdt_col = Simp_Dyn_new(t,x)
    global b_damp x_nom u_nom K_col alpha t0 uff_col nState mControl u_temp
    
    xnom(1) = interp1(t0,x_nom(1,:)',t);
    xnom(2) = interp1(t0,x_nom(2,:)',t);
   
    unom = interp1(t0,u_nom,t);

    uff_ = interp1(t0,uff_col,t);
    K_ = interp1(t0,K_col',t);
    
    xbar = x-xnom';
    ustar = unom + alpha*reshape(uff_,[mControl,1]) + reshape(K_,[mControl, nState])*(xbar);
    u_temp = [u_temp; ustar];
    dxdt = [x(2);
                ustar-b_damp*x(2)-sin(x(1))];
    dxdt_col = reshape(dxdt,[nState,1]); 
end

function [Alin, Blin] = linearize()
    global x_nom N nState mControl b_damp t0 dt
    
    Alin = zeros(nState,nState,N);
    Blin = zeros(nState,mControl,N);
%     for i = 1:N
%         Alin(:,:,i) = [0 1;-cos(x_nom(1,i)) -b_damp];
%         Blin(:,:,i) = [0;1];
%     end

    for i = 1:N
        A = [0 1;-cos(x_nom(1,i)) -b_damp];
        B = [0;1];
        M = [A B; zeros(1,nState+mControl)].*dt;
        MM = expm(M);
        Alin(:,:,i) = MM(1:nState,1:nState);
        Blin(:,:,i) = MM(1:nState,nState+1:end);
    end
end

function [uff_out, K_out] = backward_ricatti_forS(Alin, Blin)
    global q0 q1 Q P r R q0f q1f Qf S2 S1 S0 
    global x_nom u_nom x_d u_d nState mControl N
    
    End_point_weight = 2000;
    Running_weight_state = 0.5;
    Running_weight_control = 0.5;
    
    
    
    A = Alin;
    B = Blin;
    delta_x = x_d - x_nom;
    delta_u = u_d - u_nom;
    Q_ = Running_weight_state*4*Q;
    R_ = Running_weight_control*4*R;
    
    S2=zeros(nState,nState,N);
    S1=zeros(nState,N);
    S0=zeros(1,N);
    

    
    S2(:,:, N) = End_point_weight*4*Qf;
    S1(:,N) = -End_point_weight*2*Qf*delta_x(:,N);
    S0(N) = delta_x(:,N)'*Qf*delta_x(:,N);
    
    uff_out = zeros(mControl,N);
    K_out = zeros(mControl, nState, N);
    
    for i = N-1:-1:1
        q1 = -2*Q*delta_x(:,i);
        r = -2*R*delta_u(:,i);
        H = R_ + B(:,:,i)'*S2(:,:, i+1)*B(:,:,i);
        G = P + B(:,:,i)'*S2(:,:, i+1)*A(:,:,i);
        g = r + B(:,:,i)'*S1(:, i+1);% - R*delta_u(:,i) - P*delta_x(:,i);

        K = -inv(H)*G;
        uff_bar = -inv(H)*g;
        
        S2(:,:,i) = Q_ + A(:,:,i)'*S2(:,:, i+1)*A(:,:,i) + K'*H*K + 2*K'*G;
        S1(:,i) = q1 + A(:,:,i)'*S1(:, i+1) + K'*H*uff_bar + K'*g + G'*uff_bar;% - Q*delta_x(:,i) - P'*delta_u(:,i);
        
        K_out(:,:,i) = K;
        uff_out(:,i) = uff_bar;

    end
end