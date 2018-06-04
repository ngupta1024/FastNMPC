function [xdot, xNext]=augPendDynamics(x,u,modelParams)
    xdot(1,:)=x(2,:);
    xdot(2,:)=-(modelParams.g/modelParams.length)*sin(x(1,:))+...
        u/(modelParams.m*modelParams.length^2)-(((x(1,:).^2)-1).*x(2,:));
    xNext=x+xdot*modelParams.dt;
end
