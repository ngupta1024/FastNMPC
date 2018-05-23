function J=computeActualCost(actual_traj. desired_traj)
    x_diff=actual_traj.x-desired_traj.x;
    u_diff=actual_traj.u-desired_traj.u;
    J=x_diff(:,end)'*modelParams.Qf*x_diff(:,end)
end