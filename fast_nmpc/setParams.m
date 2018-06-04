function modelParams=setParams()
%% algo 1

    modelParams.g=1;%9.8
    modelParams.m=1;
    modelParams.length=1;
    modelParams.c=0;
    
    modelParams.dt=0.1;
    modelParams.T=10; %N=T/dt
    modelParams.N=modelParams.T/modelParams.dt+1;
    
    modelParams.Qt=diag([1,1]);
    modelParams.Qf=diag([10,10]);
    modelParams.Rt=1;
    
    modelParams.x_init=[0;0];
    modelParams.u_lim=1;
    
    modelParams.gen_traj=0;
    modelParams.viz=1;
    modelParams.ls_steps=1;
    modelParams.alpha_d=1.3;
    
%% algo 2

    modelParams.policy_lag=10;
    
end