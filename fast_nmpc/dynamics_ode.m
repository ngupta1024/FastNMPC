function xdot=dynamics(t,x,xO,uO,xD,uD,timestamp,K1,K2,tS1,tS2)
        xOInterp(1)=interp1(timestamp',xO(:,1),t);
        xOInterp(2)=interp1(timestamp',xO(:,2),t);
        xDInterp(1)=interp1(timestamp',xD(:,1),t);
        xDInterp(2)=interp1(timestamp',xD(:,2),t);
        uOInterp=interp1(timestamp',uO,t);
        uDInterp=interp1(timestamp',uD,t);
        K1Interp(1)=interp1(tS2,K1(:,1),t);
        K1Interp(2)=interp1(tS2,K1(:,2),t);
        K2Interp=interp1(tS1,K2',t);
        xBar=x-xOInterp';
        uDBar=uDInterp-uOInterp;
        u_=uDInterp-K1Interp*xBar-K2Interp;
        xdot=[x(2); u_ - sin(x(1))-0.3*x(2)];
    end