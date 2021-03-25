
include("read_data.jl")
include("models.jl")


hydrosys_folder = "p1"

hydrosys_data = read_data(hydrosys_folder)

model = create_model(hydrosys_data, T=24)

# T = 24			# numero total de intervalos de tempo (1 por hora)
# R = 4 			# numero de reservatorios a serem acoplados
# J = [3 3 3 5]		# numero de unidades (turbinas) para cada usina (reservatorio)	
# L = [3 5 4 2]		# numero de retas aproximando a funcao alpha (custo da agua) para cada reservatorio	
# x0 = 2000
# a = [243 1.07 -0.011 -0.000000521 -0.00000000000924]
# a_bar0 = 3*a[5]*x0^4 + a[4]*x0^3 + a[1]
# a_bar1 = - 8*a[5]*x0^3 - 3*a[4]*x0^2 + a[2]
# a_bar2 = 6*a[5]*x0^2 + 3*a[4]*x0 + a[3]
# ###############################################################################
# ############################ Parametros #######################################
# ###############################################################################
# y = []
# c1 = 
# G = 9.8066e-3
# q_min = 
# q_max = 
# v_min = [800 8000 8000 8000]			# volume minimo de cada reservatorio
# v_max = [100 100 100 100]			# volume maximo de cada reservatorio	
# R_up = [[], [] , [1,2], [3]] 			# reservatorios a montante de casa usina
# tau = []					# tempo de viagem da agua entre reservatorios
# fcm_final_dia = 0.3				# fcm percentual final (esperado) de cada reservatório
# v_inicial = 1.3*v_min				# volume inicial de cada reservatório
# preco_hora = rand(24)				# valor da energia por hora do dia
# u_begin = []					# estado inicial de cada uma das turbinas de cada reservatorio
# max_epsilon = 5					# número maimo de ligacoes/desligamento de turbinas durante o planejamento
# delta_tempo = 8					# tempo mínimo até religar a unidade (turbina)  
# pg_rampa = 10					# delta de potência (absoluta) entre 2 intervalos de tempos para cada turbina
# K = []						# número de zona da curva colina para cada turbina de cada usina
# alpha_demanda = 0.95				# multiplicador da demanda a ser cumprida
# delta_fcm = .85
# fcm_max = sum(a[i]*v_max[1]^(i-1) for i in eachindex(a)) # poly(a,v_max) 
# K_pusina = 0