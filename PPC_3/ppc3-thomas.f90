! PPC 3 - Método de Thomas
! Aluno: Igor de Oliveira Rossoni	Matrícula: 222031279
! Este programa resolve o sistema linear [A]{T} = {B} para o problema de temperaturas em um reator utilizando o método de Thomas.
! Para os cálculos sem geração interna, deve-se alterar o valor de q_dot para 0.0.
! Para calcular com geração interna, o valor utilizado para q_dot é 2.0e7.
! Os resultados do sistema linear em função do tempo estão disponíveis nos gráficos entitulados de: 
! - "campo_com_geracao.png" e "evolucao_com_geracao.png" para o caso do reator com geração interna;
! - "campo_sem_geracao.png" e "evolucao_sem_geracao.png" para o caso do reator com geração interna.

! Ao rodar este programa, espera-se um arquivo "saida.dat" que conterá as seguintes colunas: tempo (s), nó, posição x (m), temperatura (°C), respectivamente.
! Os gráficos de campo de cores foram gerados para mostrar a transferência de calor na parede do reator e os gráficos de perfil, para mostrar como a temperatura evolui em função do tempo para cada nó individualmente.


program thomas
implicit none

! Variáveis físicas
real(8) L, alpha, k_term	! L = espessura, alpha = difusividade do material, k_term = condutividade do material, 
real(8) h_conv, T_inf	! h_conv = coef. de transf. de calor por convecção, T_inf = temperatura do fluido externo
real(8) rho, Cp		! rho = densidade do material, Cp = calor específico a pressão constante
real(8) q_dot, A_fonte	! q_dot = taxa de geração, A_fonte = termo fonte agrupado

! Parâmetros da malha e tempo
integer, parameter :: n = 7 ! n = número de nós
real(8) dx, dt, t_final	! dx = espaço entre nós, dt = passo de tempo, t_final = tempo total de simulação
real(8) Fo, Bi		! Fo = número de Fourier, Bi = número de Biot
real(8) x_pos		! Posição na parede com relação ao número de nós

! Algoritmo de Thomas
real(8), dimension(n) :: e, f, g, r ! e = vetor diagonal inferior, f = vetor diagonal principal, g = vetor diagonal superior, r = vetor de fontes
real(8), dimension(n) :: e_orig, f_orig, g_orig, r_orig	 ! vetores originais
real(8), dimension(n) :: T, T_novo	 ! T = vetor de temperatura no instante atual, T_novo = vetor de temperatura no próximo passo de tempo

! Contadores
integer i, k, p, p_total   ! p = passo de tempo atual, p_total = número total de passos de tempo

open(1, file='saida.dat', status='replace', action='write')
! Considerando Dióxido de Urânio (UO2) como combustível e água pressurizada para fluido refrigerante externo
L = 0.012		! Espessura [m]
k_term = 3.0		! Condutividade térmica do UO2 [W/m-K]
rho = 10970.0		! Densidade do UO2 [kg/m3]
Cp = 300.0		! Calor específico do UO2 [J/kg-K]
alpha = k_term/(rho*Cp)	! Difusividade [m2/s]
h_conv = 5000.0		! Coef. de convecção forçada da água [W/m2-K]
T_inf = 280.0		! Temperatura do fluido refrigerante externo [°C]
!q_dot = 0.0		! Caso sem geração interna
q_dot = 2.0e7		! Caso com geração interna
dx = L/(n-1)		! Espaço entre cada nó [m]
dt = 0.1		! Passo de tempo de 0.1 s
t_final = 10.0*60.0/3.0	! Tempo final de simulação em segundos
p_total = t_final/dt	! Quantidade de passos para cada nó
Fo = alpha*dt/(dx**2)	! Cálculo do número de Fourier
Bi = h_conv*L/k_term	! Cálculo do número de Biot
A_fonte = (q_dot/(rho*Cp))*dt

! Vetor f -> diagonal principal
f_orig(1) = 1.0 + 2.0*Fo
do i = 2, n - 1
	f_orig(i) = 1.0 + 2.0*Fo
end do
f_orig(n) = 1.0 + 2.0*Fo + 2.0*Fo*Bi

! Vetor e -> diagonal inferior
e_orig(1) = 0.0
do i = 2, n - 1
	e_orig(i) = -Fo
end do
e_orig(n) = -2.0*Fo

! Vetor g -> diagonal superior
g_orig(1) = -2.0*Fo
do i = 2, n - 1
	g_orig(i) = -Fo
end do
g_orig(n) = 0.0

T = 20.0	! Temperatura inicial de 20 °C
do p = 1, p_total
	f = f_orig
	e = e_orig
	g = g_orig
	
	! Vetor r
	r(1) = T(1) + A_fonte
	do i = 2, n - 1
		r(i) = T(i) + A_fonte
	end do
	r(n) = T(n) + Fo*((2.0*Bi*T_inf) + (q_dot*(dx**2)/k_term))
	
	! Algoritmo de Thomas
	! Etapa de decomposição
	do k = 2, n
		e(k) = e(k)/f(k-1)
		f(k) = f(k) - e(k)*g(k-1)
	end do
	! Etapa de substituição progressiva
	do k = 2, n
		r(k) = r(k) - e(k)*r(k-1)
	end do
	! Etapa de substituição regressiva
	T_novo(n) = r(n)/f(n)
	do k = n-1, 1, -1
		T_novo(k) = (r(k) - g(k)*T_novo(k+1))/f(k)
	end do
	
	T = T_novo
	
	do i = 1, n
		x_pos = (i - 1)*dx
		write(1, '(F10.2 ,4X, I3, 4X, F10.5, 4X, F12.4)') p*dt, i, x_pos, T(i)
	end do
	write(1,*) 
end do
close(1)
end program thomas
