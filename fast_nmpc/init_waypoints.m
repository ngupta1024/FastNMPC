function wp=init_waypoints(modelParams)
%output wp - a vector of structures
for wp_iter=1:modelParams.num_wp
    x.state     = modelParams.states(:,wp_iter);
    x.rho_p     = modelParams.rho_p(wp_iter);
    x.t_p       = modelParams.t_p(wp_iter);
    x.weight_p    = modelParams.weight_p;
    wp{wp_iter} = x;
end
end