%% 应用Benders分解算法解决混合整数规划问题，导入的模型为Gurobi基本model
% 用来解决如下形式的规划问题:
%       BendersDecompostion (used in this function)                ( GurobiType
%         min  d'*y+c'*x                                            min   model.obj'*x + x'*model.Q*x + alpha
%         s.t. F*y + E*x <= h;
%              A*y <= b                                            s.t.  model.A*x =or<=or>=(model.sense) b (model.rhs),
%              lx<=x<=ux                                                 model.lb <= x <= model.ub,  some xj integral model.vtype)
%              ly<=y<=uy                                                 x'*Qc*x + q'*x <= beta,
%      y: binaries（complicating variables）,x(>=0): continious
%      MP1:
%         min_(z,y) z
%         s.t. z>=d'*y
%              A*y<=b
%              y: binaryies
%      SP1 (1-y\hat):
%         min_(x)  c'*x
%         s.t. E*x <= h-F*y\hat
%               x>=0
%      SP2 (dual of SP1):
%         max_(u)  (F*y\hat)'*x
%         s.t. E'*u <= c;
%               u<=0

t_BD_s = tic;  
%%  Benders decompostion parameter setting and initialize
epsilon = 1e-5;
z_UB=inf;     % Opt value of current best feasible solution
z_LB=-inf;    % Opt value of relaxed problem (MP)
MipDisplayInterval=10000;
MaxIter=100;  % Maximum iteration number
Outputflag=0;  % 0 means no MIP solving output

%% variables transform into form in BendersDecomposition variables
Ind_y=[]; % Index of y, the complicating variables, or binary variables
Ind_x=[]; % Index of x
for i=1:length(model.vtype)
    if isequal(model.vtype(i),'B')
        Ind_y=[Ind_y;i];
    end
    if isequal(model.vtype(i),'C')
        Ind_x=[Ind_x;i];
    end
end
N_y=length(Ind_y);
N_x=length(Ind_x);
N_p=N_y+N_x;    % number of variables in primal problem
N_d=length(model.sense);   % number of variables in dual problem

% In Gurobi model,
% 1) by default, model.sense are all '<', thus all coefficent matrices (F,A,E) and rhs(b,h) might need
% to be inversed
% 2) binary variables emerged before continous variables
FA=model.A(:,1:N_y);    % Both F and A are coefficients of y (compliacting variables)
E=model.A(:,(N_y+1):end);
hb=model.rhs;
d=model.obj(Ind_y);
c=model.obj(Ind_x);
MP_rows=[];   % A related rows
SP_rows=[];
for i=1:N_d
    if isempty(find(E(i,:)))
        MP_rows=[MP_rows;i];
    else
        SP_rows=[SP_rows;i];
    end
end
A=FA(MP_rows,:);
F=FA(SP_rows,:);
b=hb(MP_rows);
h=hb(SP_rows);

%% Step 1: Solve MP1 (3.4)
% if E contains some lines with all 0, then there are some constraints in
% MP too.
MP.obj=[1;zeros(N_y,1)];   % Variables in MP1 : [z;y];
MP.A=sparse([1,-d']);   %[z>=d'*y，A*y>=b （不含x的相关约束条件）]
MP.A=[MP.A;zeros(length(MP_rows),1),A];
MP.sense=['>';model.sense(MP_rows)];
MP.rhs=[0;b];
MP.lb=[0;model.lb(Ind_y)];
MP.ub=[Inf;model.ub(Ind_y)];
MP.vtype=['C';model.vtype(Ind_y)];
MP.params.DisplayInterval=MipDisplayInterval;
MP.params.outputflag=Outputflag;
r_MP=gurobi(MP,MP.params);
z_LB=r_MP.objval;

p = 0;
q = 0;

%% Construct SP1 (3.5) (1-y\hat) or SP2 (3.6) (dual of SP1)
SP.A=E(SP_rows,:);
SP.obj=c;
SP.sense=model.sense(SP_rows);
SP.lb=model.lb(Ind_x);
SP.ub=model.ub(Ind_x);
SP.vtype=model.vtype(Ind_x);
SP.params.DisplayInterval=MipDisplayInterval;
SP.params.outputflag=Outputflag;
r_SP=gurobi(SP,SP.params);
assign(recover(r_model.used_variables(Ind_x)),r_SP.x);

timeS = 0;
iter=1;
s_P_BD=[];
s_u_BD=[];
s_v_BD=[];
s_w_BD=[];
abs_error=abs((z_UB-z_LB)/z_UB);

while iter<=MaxIter
    %% Step 2: Solve the SP1 (3.5) (1-y\hat) or SP2 (3.6) (dual of SP1)
    y_hat=r_MP.x(2:end);
    SP.rhs=h-F*y_hat;
    r_SP=gurobi(SP,SP.params);

    if strcmp(r_SP.status,'OPTIMAL')
        p = p+1;
        display(['Add optimality cut ',num2str(p),' !']);
        assign(recover(r_model.used_variables(Ind_x)),r_SP.x);
        s_P_BD((3*p-2):3*p,:)=value(P);
        %% add optimality cut
        z_UB=r_SP.objval+d'*r_MP.x(2:end);
        abs_error=abs((z_UB-z_LB)/z_UB);
        if abs_error<=epsilon
            display(['Upper Bound: ', num2str(z_UB),'  Lower Bound: ', num2str(z_LB),'  Gap: ',num2str(round(abs_error*100,2)),'%']);
            break
        end
        pi=r_SP.pi;
        oc(p).A = [1,(pi'*F-d')];  % optimality cut: z>=d'*y+(h-F*y)'*pi → z+(pi'*F-d')*y>=pi'*h
        oc(p).b = pi'*h;
        MP.A=[MP.A;oc(p).A];
        MP.rhs=[MP.rhs;oc(p).b];
        MP.sense=[MP.sense;'>'];
    elseif strcmp(r_SP.status,'INF_OR_UNBD')
        % Calculate u^r in SP1(3.7)
        SP_n=SP;
        SP_n.obj=[zeros(N_x,1);ones(length(SP_rows),1)];
        SP_n.A=[SP.A,-diag(ones(length(SP_rows),1))];
        SP_n.lb=[SP.lb;zeros(length(SP_rows),1)];
        SP_n.ub=[SP.ub;ones(length(SP_rows),1)*Inf];
        SP_n.vtype=[SP.vtype;ones(length(SP_rows),1)*'C'];
        r_SP_n=gurobi(SP_n,SP.params);
        q = q+1;
        display(['Add feasibility cut ',num2str(q),' !']);
        pi=r_SP_n.pi;
        fc(q).A = [0,pi'*F];  % feasibility cut: (h-F*y)'*u<=0 → u'*F*y>=u'*h
        fc(q).b = pi'*h;
        MP.A=[MP.A;fc(q).A];
        MP.rhs=[MP.rhs;fc(q).b];
        MP.sense=[MP.sense;'>'];
    end
    %% Step 3: Solve MP2 to obtain a new lower bound solution z_LB w.r.t. y_hat
    r_MP=gurobi(MP,MP.params);
    assign(recover(r_model.used_variables(Ind_y)),r_MP.x(2:end)); % exclude varialbe z in MP
    s_u_BD((3*iter-2):3*iter,:)=value(u);
    s_v_BD((3*iter-2):3*iter,:)=value(v);
    s_w_BD((3*iter-2):3*iter,:)=value(w);
    z_LB=r_MP.objval;
    iter=iter+1;
    abs_error=abs((z_UB-z_LB)/z_UB);
    display(['Upper Bound: ', num2str(z_UB),'  Lower Bound: ', num2str(z_LB),'  Gap: ',num2str(round(abs_error*100,2)),'%']);
end
t_BD_e = toc(t_BD_s);  
display(['采用Gurobi+benders分解所用计算时间: ',num2str(round(t_BD_e,2)),' s']);

