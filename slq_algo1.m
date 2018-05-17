function slq_algo1
modelParams=setParams();


%% Given
    function modelParams=setParams()
        modelParams.g=9.8;
        modelParams.length=1;
        modelParams.dt=0.01;
    end

    function xNext=simplePend_dynamics(x,u, modelParams)
        xdot(1)=x(2);
        xdot(2)=-(modelParams.g/modelParams.length)*sin(x(1));
        xNext=x+xdot*modelParams.dt;
    end

end