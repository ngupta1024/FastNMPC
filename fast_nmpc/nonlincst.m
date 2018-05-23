function [c,ceq]=nonlincst(nom_traj,dynamics,x_goal,modelParams)
    xNext=nom_traj(2:3,:);
    u=nom_traj(1,:);
    BoundaryFinal=x_goal-xNext(:,end);
    BoundaryInit=modelParams.x_init-xNext(:,1);
    delta=hermiteSimpsonDefects(xNext,u,dynamics,modelParams);
    c = [];  %No inequality constraints
    ceq = [BoundaryInit; reshape(delta,2*modelParams.N-2,1); BoundaryFinal; u(end)];
end