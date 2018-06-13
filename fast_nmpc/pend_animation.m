function pend_animation(x,fig_pend)
figure(3);
hold on;
modelParams=setParams();
O=[0 0];
axis(gca,'equal');
axis([-modelParams.length*1.3 modelParams.length*1.3 -modelParams.length*1.3 modelParams.length*1.3]);
grid on;

% Loop for animation

for i=1:length(x)
    %Mass Point
    P=modelParams.length*[sin(x(i)) -cos(x(i))];
    
    %circle in origin
    O_circ=viscircles(O,0.01,'color','black');
    
    % pendulum
    pend=line([O(1) P(1)],[O(2) P(2)],'LineWidth',5);
    
    %ball
    ball=viscircles(P,0.03);
    
    drawnow
    
    if i<length(x)
        delete(pend);
        delete(ball);
        delete(O_circ);
    end
end
hold off;
end