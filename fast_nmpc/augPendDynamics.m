function [xdot, xNext]=augPendDynamics(x,u,modelParams)
    xdot=[x(2,:);-sin(x(1,:))+u-(((x(1,:).^2)-1).*x(2,:))];
    xNext=x+xdot*modelParams.dt;
end
