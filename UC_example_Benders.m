% UC example, data most from [1]
% [1] https://yalmip.github.io/example/unitcommitment/
% [2] G. Morales-Espa�0�9a, J. M. Latorre and A. Ramos, "Tight and Compact MILP Formulation for the Thermal Unit Commitment Problem," in IEEE Transactions on Power Systems, vol. 28, no. 4, pp. 4897-4908, Nov. 2013, doi: 10.1109/TPWRS.2013.22514
clear all;
% close all;
clc

%% Basic condition
Nunits = 3;
Hours = 48;

Pmax = [100;50;25];
Pmin = [20;40;1];

C = [10 20 20];
Cu_NL = [2 3 4]; % No-load cost of unit
Cup = [5 10 5]; % Startup cost of unit
Cdown = [5 10 5]; % Shutdown cost of unit

Pforecast = 100 + 50*sin((1:Hours)*2*pi/24);

%% Variable definition
u = binvar(Nunits,Hours,'full');
v = binvar(Nunits,Hours,'full');
w = binvar(Nunits,Hours,'full');
P = sdpvar(Nunits,Hours,'full');

%% Constraints
Cons = [];
for t = 1:Hours
    Cons = [Cons, u(:,t).*Pmin <= P(:,t) <= u(:,t).*Pmax];
end

%% Adding minimum up- and down-time
TUg = [6;30;1];  %from minup in [1]
TDg = [3;6;3];   %from mindown in [1]

%% Cons (6) and (7) in [2]
for g=1:3
    for t=TUg(g):Hours
        Cons = [Cons,sum(v(g,(t-TUg(g)+1):t),2)<=u(g,t)];
    end
    for t=TDg(g):Hours
        Cons = [Cons,sum(w(g,(t-TDg(g)+1):t),2)<=1-u(g,t)];
    end
end

%% Cons (8) in [2]
for t=2:Hours
    Cons = [Cons,u(:,t)-u(:,t-1)==v(:,t)-w(:,t)];
end

for t = 1:Hours
    Cons = [Cons, sum(P(:,t)) >= Pforecast(t)];
    Cons = [Cons, sum(u(:,t).*Pmax) >= Pforecast(t)]; % !!! Important: Redundant constraints to limit the "complicating variables" y
end

%% Objective function
Obj_f = 0;
Obj_u = 0;
Obj_up = 0;
Obj_down = 0;
for t = 1:Hours
    Obj_f=Obj_f+C*P(:,t);
    Obj_u=Obj_u+Cu_NL*u(:,t);
    Obj_up=Obj_up+Cup*v(:,t);
    Obj_down=Obj_down+Cdown*w(:,t);
end
Obj=Obj_f+Obj_u+Obj_up+Obj_down;

%% (1) Solve the UC problem by Gurobi
ops = sdpsettings('solver','gurobi','verbose',1);
result_NoBD=optimize(Cons,Obj,ops);
s_P=value(P);
s_u=value(u);
s_v=value(v);
s_w=value(w);
s_Obj=value(Obj);

figure
h1=bar(s_P','stack');
legend('Unit 1','Unit 2','Unit 3');
title('UC result solved by Gurobi');

display(['Solved by gurobi，without Benders Decomposition: ', num2str(result_NoBD.solvertime),' s']);

%% (2) Solve the Problem by Benders Decomposition
% 2.1 Export Model with Gurobi
ops = sdpsettings('solver','gurobi','verbose',1);
[model,r_model,diagnostic,internalmodel] = export(Cons,Obj,ops);

% 2.2 Solve it by my BendersDecomposition (embedded with Gurobi)
BendersDecomposition_Gurobi_2022;
