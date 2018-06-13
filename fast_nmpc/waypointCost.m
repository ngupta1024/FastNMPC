function W=waypointCost(x,t,wp,modelParams)
W=0;
for wp_iter=1:modelParams.num_wp
    xcap=x-wp{wp_iter}.state;
    xcap(1,:)=wrapToPi(xcap(1,:));
    W=W+xcap'*wp{wp_iter}.weight_p*xcap*sqrt(wp{wp_iter}.rho_p/(2*pi))*exp(-wp{wp_iter}.rho_p/2*(t-wp{wp_iter}.t_p)^2);
end
end