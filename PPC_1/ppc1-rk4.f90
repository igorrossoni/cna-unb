! PPC #1 - Implementação computacional do método Runge-Kutta de quarta ordem para solução do problema de sedimentação de uma esfera em regime de baixo Reynolds.

! Aluno: Igor de Oliveira Rossoni
! Matrícula: 222031279

! Para o caso de Re = 0, observa-se que a solução numérica obtida pelo método RK4 coincide praticamente com a solução analítica ao longo de todo o intervalo de tempo analisado. As diferenças entre os valores aparece a partir da sexta casa decimal, sendo desprezível do ponto de vista prático. Graficamente, as curvas numérica e analítica se sobrepõem completamente.

! O refinamento temporal, com a redução do passo de tempo, melhora significativamente a qualidade da solução numérica. Fazendo teste para h = 1, encontra-se um erro mais elevado, pois o passo é relativamente grande, reduzindo a capacidade do método de capturar a variação de velocidade ao longo do tempo, desse modo, a solução pode apresentar desvios em relação ao comportamento esperado. Ao reduzir o passo para h = 0.5, é perceptível a diferença na precisão, diminuindo o erro e, graficamente, tornando a curva mais suave. Com h = 0.1, o refinamento temporal é suficiente para garantir alta precisão, fazendo com que a solução numérica se aproxime muito da solução analítica, reduzindo o erro consideravelmente.

! Para Re != 0, não há solução numérica simples em função do tempo, mas é possível determinar analiticamente a velocidade terminal impondo regime permanente, usando:
! v_analitica = (-1.0 + sqrt(1.0 + (3.0/2.0)*Res)) / ((3.0/4.0)*Res)
! A solução numérica obtida pelo método apresenta comportamento físico consistente, para Re = 0, a solução tende a 1, enquanto para Re != 0, converge para a velocidade terminal, que é menor que 1. Isso mostra confiabilidade no método RK4 e mostrando que ele consegue representar corretamente tanto o comportamento transitório quanto o regime permanente do sistema.

! Para Re = 0, o sistema tende ao regime de Stokes, em que a solução é linear e a velocidade terminal é igual a 1. Conforme o número de Reynols aumenta, surge o termo não linear associado ao arrasto quadrático, diminuindo a velocidade terminal e afastando do caso linear. Assim, quanto maior o valor de Re, maior é o desvio em relação ao regime assintótico Re = 0.

program rk4_sedimentacao
implicit none

real h, t0, tf  ! Passo de tempo, Tempo inicial, Tempo final
real t, v       ! Variáveis adimensionais de tempo e velocidade
real k1, k2, k3, k4 ! Variáveis para o método RK4
real St, Res    ! Números de Stokes e Reynolds
integer i, n
real v_analitica, erro

! Parâmetros do problema
St = 1.0    ! Número de Stokes
Res = 1.0   ! Número de Reynolds
t0 = 0.0    ! Tempo inicial
tf = 10.0   ! Tempo final
h = 0.1     ! Passo de tempo
v = 0.0     ! Condição inicial v*(0) = 0
t = t0      ! Tempo inicial é t0
n = int((tf-t0)/h)  ! Número de iterações

open(1,file='saida.dat', status='replace')  ! Abrindo o arquivo para escrever os dados resultantes dos cálculos

! Loop RK4
do i=0, n-1

k1 = f(v, St, Res)
k2 = f(v + h*k1/2.0, St, Res)
k3 = f(v + h*k2/2.0, St, Res)
k4 = f(v + h*k3, St, Res)

v = v + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)   ! Calculando o novo valor de v
t = t + h   ! Calculando o novo valor de t

if (Res == 0.0) then
v_analitica = 1.0 - exp(-t/St)  ! Solução analítica para Res = 0
erro = abs(v - v_analitica)
write(1,*) t, v, v_analitica, erro

else
v_analitica = (-1.0 + sqrt(1.0 + (3.0/2.0)*Res)) / ((3.0/4.0)*Res) ! Velocidade terminal analítica para Res != 0
erro = abs(v - v_analitica)
write(1,*) t, v, v_analitica, erro
end if

end do

! Definindo a função
contains
real function f(v, St, Res)
implicit none
real, intent(in) :: v, St, Res
f = (1.0/St)*(1.0 - v - (3.0/8.0)*Res*v*v)  ! Equação adimensionalizada da sedimentação de uma partícula
end function f

end program rk4_sedimentacao
