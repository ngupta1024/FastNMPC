%% dynamics of a simple pendulum
% params x: state vector 2x1
% params u: input vector 1
% params modelParams: struct
%returns xNext: 2x1, xdot=2x1
function [x_dot, x_next]=simplePendDynamics(x,u, modelParams)
    x_dot(1,:)=x(2,:);
    x_dot(2,:)=-(modelParams.g/modelParams.length)*sin(x(1,:))-((modelParams.c/modelParams.m)*x(2,:)) +(u/(modelParams.m*modelParams.length^2));    
    x_next=x+x_dot*modelParams.dt;
end