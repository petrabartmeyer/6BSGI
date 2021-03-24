using JuMP, Gurobi, DataFrames

T = 24			# numero total de intervalos de tempo (1 por hora)
R = 4 			# numero de reservatorios a serem acoplados
J = [3 3 3 5]		# numero de unidades (turbinas) para cada usina (reservatorio)	
L = [3 5 4 2]		# numero de retas aproximando a funcao alpha (custo da agua) para cada reservatorio	
a_bar0 = 3*a4*x0^4 + a3*x0^3 + a0;
a_bar1 = - 8*a4*x0^3 - 3*a3*x0^2 + a1;
a_bar2 = 6*a4*x0^2 + 3*a3*x0 + a2;
###############################################################################
############################ Parametros #######################################
###############################################################################
y = []
c1 = 
G = 9.8066e-3
q_min = 
q_max = 
v_min = [800 8000 8000 8000]			# volume minimo de cada reservatorio
v_max = [100 100 100 100]			# volume maximo de cada reservatorio	
R_up = [[], [] , [1,2], [3]] 			# reservatorios a montante de casa usina
tau = []					# tempo de viagem da agua entre reservatorios
fcm_final_dia = 0.3				# fcm percentual final (esperado) de cada reservatório
v_inicial = 1.3*v_min				# volume inicial de cada reservatório
preco_hora = rand(24)				# valor da energia por hora do dia
u_begin = []					# estado inicial de cada uma das turbinas de cada reservatorio
max_epsilon = 5					# número maimo de ligacoes/desligamento de turbinas durante o planejamento
delta_tempo = 8					# tempo mínimo até religar a unidade (turbina)  
pg_rampa = 10					# delta de potência (absoluta) entre 2 intervalos de tempos para cada turbina
K = []						# número de zona da curva colina para cada turbina de cada usina
alpha_demanda = 0.95				# multiplicador da demanda a ser cumprida
delta_fcm = .85
fcm_max =  poly(a,v_max) 
K_pusina = 0

function modelo(data::Union{DataFrame, Dict})
    model = Model(Gurobi.Optimizer)
    ###############################################################################
    ############################ Funcoes auxiliares ###############################
    ###############################################################################
    poly(a::Vector,x::Float64) = sum(a[i]*x^(i-1) for i in eachindex(a)) 
    pmt(pg,g) =  poly(g,pg)
    pgg(pg,f) = f[1]*exp(f[2]*pg)
    h(a,b,k_p,k_pusina,k_s,q,V,s) = poly(a,V) - poly(b,Q+s) + k_p*q^2 + K_pusina*Q^2 + k_s*q^2
    rho(c,q,h) = c[1] + c[2]*q +c[3]*h + c[4]*h*q + c[5]*q^2+ c[6]*h^2 
    pst(G,c,q,h) = G*rho(c,q,h)*h*q


    ###############################################################################
    ############################ Variáveis ########################################
    ###############################################################################
    @variable(model, 0<= v[r=1:R, t=1:T])			# volume armazenado no reservatorio R no tempo t
    @variable(model, 0<= s[r=1:R, t=1:T])			# volume vertido no reservatorio R no tempo t
    @variable(model, 0<= Q[r=1:R, t=1:T])			# vazao turbinada no reservatorio R no tempo t
    @variable(model, 0<= q[j=1:max(J), r=1:R, t=1:T])		# vazao turbinada pela unidade j do reservatorio R no tempo t
    @variable(model, 0<= pg[j=1:max(J), r=1:R, t=1:T])		#u_jrt = 1 se a unidade j do reservatorio r no tempo T está ligado
    #@variable(model, z[j=1:max(J), k=1:phi, r=1:R, t=1:T], Bin)     # controlar curvas colina
    #@variable(model, u[j=1:max(J), r=1:R, t=1:T], Bin)		#u_jrt = 1 se a unidade j do reservatorio r no tempo T está ligado
    #@variable(model, u_on[j=1:max(J), r=1:R, t=1:T], Bin)	#u_jrt = 1 se a unidade j do reservatorio r foi ligada no tempo T  
    #@variable(model, u_off[j=1:max(J), r=1:R, t=1:T], Bin)	#u_jrt = 1 se a unidade j do reservatorio r  foi desligada no tempo T
    
    ###############################################################################
    ############################ Objective function################################
    ###############################################################################
  
    @objective(model, sum(preco_hora[t]*pg[r,t] for r=1:R, t=1:T)) # + sum(alpha[r] for r=1:R))
    
    ###############################################################################
    ############################ Variáveis ########################################
    ###############################################################################
    @constraint(model, [r=1:R, t=1], v[r,t] - v_inicial[r] + c1*(Q[r,t]+s[r,t]) == c1*y[r,t] )
    @constraint(model, [r=1:R, t=2:T], v[r,t] - v[r, t-1] + c1*(Q[r,t]+s[r,t] - sum(Q[m,t-tau[m,r]]+s[m,t-tau[m,r]] for m in R_up[r] if (t-tau[m,r])>=1) == c1*y[r,t] )
    @constraint(model, [r=1:R, t=1:T], v_min[r] <= v[r,t])
    @constraint(model, [r=1:R, t=1:T], v_max[r] >= v[r,t])
    @constraint(model, [r=1:R, t=1:T], sum(pg[j,r,t] for j=1:n[r,t]) >= alpha_demanda*L[r,t])
    @constraint(model, [r=1:R, t=1:T, j=1:J[r]], pg[j,r,t] - pst[j,r,t] + pmt[j,r,t] + pgg[j,r,t] ==0)
    @constraint(model, [r=1:R, t=1:T], sum(q[j,r,t] for j=1:J[r])-Q[r,t]==0)
    @constraint(model, [r=1:R, t=1:T, j=1:J[r]], q_min[r,j,t] <= q[j,r,t])
    @constraint(model, [r=1:R, t=1:T, j=1:J[r]], q_max[r,j,t] >= q[j,r,t])
    @constraint(model, [r=1:R, t=1:T, j=1:J[r]], pg_min[r,j,t] <= pg[j,r,t])
    @constraint(model, [r=1:R, t=1:T, j=1:J[r]], pg_max[r,j,t] >= pg[j,r,t])
   # @constraint(model, [r=1:R, t=1:T, j=1:J[r]], pg_min[r,j,t]*z[j,k,r,t] <= pg[j,r,t])
   # @constraint(model, [r=1:R, t=1:T, j=1:J[r],k=1:K], pg_max[r,j,t]*z[j,r,k,t] >= pg[j,r,t])
    
    ###############################################################################
    ############################ FCM ao final do dia ##############################
    ############################ fcm = a0+ a1*v+a2v^2 #############################
    ###############################################################################
    @constraint(model, [r=1:R, t=T], poly(a,v[r,t]) >= fcm_max[r]*delta_fcm[r])	
   # @constraint(model, [r=1:T, t=T, l=1:L[r]], alpha[r] >= coef[r,l,1] + coef[r,l,2]*v[r,t])
    
    
    ###############################################################################
    ######################### Rampas para  das turbinas ####################
    ###############################################################################
    @constraint(model, [r=1:R, t=1:T-1, j=1:J[r]], pg[j,r,t]-pg[j,r,t+1]<=pg_rampa[j,r]*pg[j,r,t])    
    @constraint(model, [r=1:R, t=1:T-1, j=1:J[r]], pg[j,r,t+1]-pg[j,r,t]<=pg_rampa[j,r]*pg[j,r,t])
    
    ###############################################################################
    ######################### ativacao/destivacao das turbinas ####################
    ###############################################################################
   # @constraint(model, [r=1:R, t=1, j=1:J[r]], u_on[j,r,t] >= u[j,r,t] - u_begin[j,r])
   # @constraint(model, [r=1:R, t=1, j=1:J[r]], u_off[j,r,t] >= u_begin[j,r] - u[j,r,t])
   # @constraint(model, [r=1:R, t=2:T, j=1:J[r]], u_on[j,r,t] >= u[j,r,t] - u[j,r,t-1])
   # @constraint(model, [r=1:R, t=2:T, j=1:J[r]], u_off[j,r,t] >= u[j,r,t-1] - u[j,r,t])
   # @constraint(model, [r=1:R, t=1:T-delta_tempo, j=1:J[r]], delta_tempo*u_off[r,j,t] <= delta_tempo -sum(u[r,j,t+delta_tempo]))
   # @constraint(model, sum(u_on[j,r,t]+u_off[j,r,t] for r=1:R, t=1, j=1:J[r]) <= max_epsilon)
    
    ###############################################################################


    ###############################################################################
    ############## restricoes se usar as curvas colina ############################
    ###############################################################################
    #@constraint(model, [r=1:R, t=1:T, j=1:J[r]], pg_min[r,j,t]*z[j,r,t] <= pg[j,r,t])
    #@constraint(model, [r=1:R, t=1:T, j=1:J[r]], pg_max[r,j,t]*z[j,r,t] >= pg[j,r,t])
    #@constraint(model, [r=1:R, t=1:T, j=1:J[r]], z[r,j,t] == u[j,r,t])
    #@constraint(model, [r=1:R, t=1:T, j=1:J[r]], z[r,j,t] <= u[j,r,t])

    return model
end
