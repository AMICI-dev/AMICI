function [model] = model_calvetti_syms()

% set the parametrisation of the problem options are 'log', 'log10' and

model.param = 'lin';

%%
% STATES
% create state syms 
syms V1 V2 V3 f1 f2 f3
% create state vector
model.sym.x = [V1, V2, V3, f1, f2, f3];
%%
% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms                                                                
% create parameter vector 
model.sym.p = [ ];
%%  
% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
syms V1ss R1ss V2ss R2ss V3ss R3ss
% create parameter vector 
model.sym.k = [V1ss, R1ss, V2ss, R2ss, V3ss, R3ss];
%%
% SYSTEM EQUATIONS
% create symbolic variable for time
syms t f0 
model.sym.xdot = sym(zeros(size(model.sym.x)));
p1=1;
p2=1-R1ss;
p3=1-(R1ss+R2ss);
L1=(R1ss*(V1ss^2))^(1/3);
C1ss=V1ss/(p1-0.5*R1ss);

C2ss=V2ss/(p2-0.5*R2ss);
L2=(R2ss*(V2ss^2))^(1/3);
C3ss=V3ss/(p3-0.5*R3ss);
L3=(R3ss*(V3ss^2))^(1/3);
R2=(L2^3)/(V2^2);
R1=(L1^3)/(V1^2);
R3=(L3^3)/(V3^2);
s=am_stepfun(t,10,1,12,0);
model.sym.xdot(1)=(1/31)*(((1.29 -(V1/V1ss))/(1.29-1))+s-2*(V1/(C1ss*((R1+R2)*f1+(R2+R3)*f2+R3*f3))));
model.sym.xdot(2)=(1/163)*(((1.51 -(V2/V2ss))/(1.51-1))- 2*(V2/(C2ss*((R2+R3)*f2+R3*f3))));
model.sym.xdot(3)=(1/122)*(((1000000 -(V3/V3ss))/(1000000-1))- 2*(V3/(C3ss*(R3*f3))));
f0=(2/R1)-((R1+R2)*f1 + (R2+R3)*f2 +R3*f3)/R1;
model.sym.xdot(4)=f0-model.sym.xdot(1)-f1;
model.sym.xdot(5)=f1-model.sym.xdot(2)-f2;
model.sym.xdot(6)=f2-model.sym.xdot(3)-f3;


model.sym.M=[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;];
%%
% INITIAL CONDITIONS
model.sym.x0(1) = V1ss;
model.sym.x0(2)=V2ss;
model.sym.x0(3)=V3ss;
model.sym.x0(4)=1;
model.sym.x0(5)=1;
model.sym.x0(6)=1;
model.sym.dx0 = sym(zeros(size(model.sym.x)));

% OBSERVALES
model.sym.y = sym(zeros(1,6));
model.sym.y(1) =V1;
model.sym.y(2)=V2;
model.sym.y(3)=V3;
model.sym.y(4)=f0;
model.sym.y(5)=f1;
model.sym.y(6)=f2;
end

