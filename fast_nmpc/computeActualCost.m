function J=computeActualCost(actual_traj,desired_traj,modelParams)
    wp=init_waypoints(modelParams);
    x_diff=actual_traj.x-desired_traj.x;
    x_diff(1,:)=wrapToPi(x_diff(1,:));
    u_diff=actual_traj.u-desired_traj.u;
    J=x_diff(:,end)'*modelParams.Qf*x_diff(:,end);
    for inst_iter=1:modelParams.N-1
        J=J+x_diff(:,inst_iter)'*modelParams.Qt*x_diff(:, inst_iter)...
            +u_diff(inst_iter)'*modelParams.Rt*u_diff(inst_iter);
        if modelParams.wp_bool
            J=J+waypointCost(actual_traj.x(:,inst_iter),inst_iter*modelParams.dt,wp,modelParams);
        end
    end    
end