%% params nomTraj=[u;x] %3xN
%% params modelParams 
%% return J scalar cost function
function J=costFunction(nomTraj, modelParams)
%used in generateTraj.m

%%only instantaneous term
J=0;
for iter=1:modelParams.N-1
    J=J+nomTraj(2:3,iter)'*modelParams.Qt*nomTraj(2:3,iter)...
        +nomTraj(1,iter)'*modelParams.Rt*nomTraj(1,iter);
end

%final cost
% J=J+nomTraj(2:3,modelParams.N)'*modelParams.Qf*nomTraj(2:3,modelParams.N);

end