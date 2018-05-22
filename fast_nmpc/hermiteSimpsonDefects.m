function delta=hermiteSimpsonDefects(xNext,u, dynamics,modelParams)
xLow = xNext(:,1:(end-1));
xHi = xNext(:,2:end);
[xdot,~] = dynamics(xNext,u);  %no time dependence
xdotLow = xdot(:,1:(end-1));
xdotHi = xdot(:,2:end);
uLow=u(:,1:end-1);
uHi=u(:,2:end);
xck = 0.5*(xLow + xHi)+ (modelParams.dt/8).*(xdotLow-xdotHi);
uck=0.5*(uLow+uHi);
% xdotck=[xck(2,:);-sin(xck(1,:))+uck-(((xck(1,:).^2)-1).*xck(2,:))];
[xdotck,~]=dynamics(xck,uck);
delta = xLow-xHi + (modelParams.dt/6)*(xdotLow + 4*xdotck + xdotHi);
end