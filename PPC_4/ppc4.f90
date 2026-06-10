! PPC 4 - Maximização de duas variáveis
! Aluno: Igor de Oliveira Rossoni       Matrícula: 222031279

! Este programa resolve um problema de maximização de duas variáveis por dois métodos diferentes, por Aclive Máximo (AM) e por Gradientes Conjugados (Fletcher-Reeves) (GC). A função a ser maximizada é f(x,y) = 2xy + 2x - x² - 2y², cujo ótimo analítico é (x, y) = (2, 1). O código foi validado com os pontos iniciais (x0, y0) = (-2, 3), mas funciona para qualquer ponto inicial de busca.

! O método do aclive máximo utiliza a direção pelo gradiente da função enquanto o método dos gradientes conjugados atualiza a direção de busca pelo gradiente da função no ponto atual somado a um parâmetro minimizador multiplicado pelo passo de busca anterior. O parâmetro minimizador é dado pela norma do gradiente atual ao quadrado divido pela norma do gradiente anterior ao quadrado.

! A saída esperada desse programa são três arquivos, "output1.dat" e "output2.dat", que são os logs dos métodos do AM e do GC respectivamente, e "function.dat", que são amostras de (x, y, f) para plotagem das curvas de nível. Após os resultados, deve-se rodar o programa "grafico.py" para gerar o gráfico com as curvas de nível e a comparação entre os dois métodos.

! Para calcular o passo dos métodos, foi utilizada interpolação quadrática na linha de busca, utilizando a fórmula do vértice, e, caso a parábola for degenerada, é adotado como passo o melhor entre h0, h1 e h2. Fazendo a análise do gráfico, percebe-se que o método AC demora muito mais para convergir do que o método GC. Além disso, percebe-se a diferença na trajetória, com o AM sendo em zigue-zague e o GC mais direto.

program ppc4
implicit none
! Variáveis
real(8) :: x, y, dx, dy, f_val 		! x, y = par ordenado | dx, dy = componentes do grad | f_val = função avaliada no par ordenado
real(8) :: passo, erro			! passo = passo de distância | erro = módulo do gradiente
real(8) :: x_inicial, y_inicial		! Par ordenado inicial de busca
real(8) :: dx_ant, dy_ant, px, py, beta	! dx_ant, dy_ant = componentes do grad no passo anterior | px, py = direções de busca | beta = parâmetro minimizador
real(8) :: x_malha, y_malha, f_malha	! Parâmetros para a geração da malha
integer :: iter
real(8), parameter :: tol = 1.0e-6

print *, "=========================================="
print *, " PPC4 - Maximização de duas variáveis"
print *, "=========================================="
print *, "Digite o valor inicial para x0:"
read *, x_inicial
print *, "Digite o valor inicial para y0:"
read *, y_inicial
print *, " "

! Método do aclive máximo
print *, "=========================================="
print *, " Metodo do Aclive Maximo"
print *, " (Steepest Accent)"
print *, "=========================================="

open(unit=10, file='output1.dat', status='replace')

passo = 0.0d0			! Passo inicial
iter = 1			! Iteração 1
x = x_inicial
y = y_inicial
f_val = funcao(x, y)		! Valor da função
dx = calc_dfx(x, y)		! Gradiente em x
dy = calc_dfy(x, y)		! Gradiente em y
erro = sqrt(dx**2 + dy**2) 	! Módulo do gradiente

write(*, '("Iter: ", I3, " | x: ", F15.8, " | y: ", F15.8, " | f: ", F15.8, " | passo (h): ", F15.8, " | erro: ", E10.3)') &
iter, x, y, f_val, passo, erro
write(10, *) iter, erro, passo, x, y, dx, dy

do while (erro.gt.tol)
        ! Interpolação quadrática
        call interpolacao(x, y, dx, dy, passo)
        
        ! Atualização dos valores de x, y e cálculo do próximo gradiente
        x = x + passo*dx
        y = y + passo*dy
        f_val = funcao(x, y)		! Valor da função
	dx = calc_dfx(x, y)		! Gradiente em x
	dy = calc_dfy(x, y)		! Gradiente em y
	erro = sqrt(dx**2 + dy**2) 	! Módulo do gradiente
	iter = iter + 1
	
	write(*, '("Iter: ", I3, " | x: ", F15.8, " | y: ", F15.8, " | f: ", F15.8, " | passo (h): ", F15.8, " | erro: ", E10.3)') &
	iter, x, y, f_val, passo, erro
	write(10, *) iter, erro, passo, x, y, dx, dy
end do

close(10)
print *, " "
write(*, *) "Log salvo em output1.dat"
print *, " "

! Método dos gradientes conjugados
print *, "=========================================="
print *, " Metodo dos Gradientes Conjugados"
print *, " (Fletcher-Reeves)"
print *, "=========================================="

open(unit=20, file='output2.dat', status='replace')

passo = 0.0d0
iter = 1
x = x_inicial
y = y_inicial
f_val = funcao(x, y)		! Valor da função
dx = calc_dfx(x, y)		! Gradiente em x
dy = calc_dfy(x, y)		! Gradiente em y
erro = sqrt(dx**2 + dy**2)	! Módulo do gradiente
px = dx				! Direção inicial de busca em x
py = dy				! Direção inicial de busca em y

write(*, '("Iter: ", I3, " | x: ", F15.8, " | y: ", F15.8, " | f: ", F15.8, " | passo (h): ", F15.8, " | erro: ", E10.3)') &
iter, x, y, f_val, passo, erro
write(20, *) iter, erro, passo, x, y, dx, dy

do while(erro.gt.tol)	
	! Interpolação quadrática
        call interpolacao(x, y, px, py, passo)
        
        ! Atualização dos valores de x e y e iteração
        x = x + passo*px
        y = y + passo*py
        f_val = funcao(x, y)		! Valor da função
        iter = iter + 1
	
	! Armazenamento e atualização dos valores do gradiente
	dx_ant = dx
        dy_ant = dy
        dx = calc_dfx(x, y)
        dy = calc_dfy(x, y)
        erro = sqrt(dx**2 + dy**2)
        
        write(*, '("Iter: ", I3, " | x: ", F15.8, " | y: ", F15.8, " | f: ", F15.8, " | passo (h): ", F15.8, " | erro: ", E10.3)') &
	iter, x, y, f_val, passo, erro
	write(20, *) iter, erro, passo, x, y, dx, dy
	
	if (erro.gt.tol) then
		beta = (dx**2 + dy**2) / (dx_ant**2 + dy_ant**2)	! Cálculo do parâmetro minimizador
		! Atualização das direções de busca
		px = dx + beta*px
		py = dy + beta*py
	end if
end do

close(20)
print *, " "
write(*, *) "Log salvo em output2.dat"
print *, " "

! Geração das amostras de x, y e f(x,y) para plotar as curvas de nível
open(unit=30, file='function.dat', status='replace')

x_malha = -2.5d0			! Começo da topografia em x
do while (x_malha.le.2.5d0)		! Fim da topografia em x
	y_malha = -1.0d0		! Começo da topografia em y
	do while (y_malha.le.3.5d0)	! Fim da topografia em y
		f_malha = funcao(x_malha, y_malha)	! Valores da função para cada ponto da malha
		write(30, *) x_malha, y_malha, f_malha
		y_malha = y_malha + 0.1d0
	end do
	write(30, *)
	x_malha = x_malha + 0.1d0
end do
close(30)
print *, " "
print *, "Arquivo function.dat gerado com sucesso"
print *, " "

contains
	real(8) function funcao(x, y) ! Função objetivo
		real(8), intent(in) :: x, y
		funcao = 2.0d0*x*y + 2.0d0*x - x**2 - 2.0d0*(y**2)
	end function funcao
	
	real(8) function calc_dfx(x, y) ! Derivada parcial com relação a x
		real(8), intent(in) :: x, y
		calc_dfx = 2.0d0*y + 2.0d0 - 2.0d0*x
	end function calc_dfx
	
	real(8) function calc_dfy(x, y) ! Derivada parcial com relação a y
		real(8), intent(in) :: x, y
		calc_dfy = 2.0d0*x - 4.0d0*y
	end function calc_dfy
	
	subroutine interpolacao(x_atual, y_atual, dir_x, dir_y, passo) ! Sub-rotina de interpolação quadrática
		implicit none
		real(8), intent(in) :: x_atual, y_atual, dir_x, dir_y	! Entrada
		real(8), intent(out) :: passo				! Saída
		real(8) :: h0, h1, h2, g0, g1, g2, num, den, g_max	! Parâmetros para interpolação
		
		! Definindo pontos iniciais para h
      		h0 = 0.0d0
        	h1 = 0.2d0
        	h2 = 0.4d0
        	! Avaliando a função ao longo da direção do gradiente
        	g0 = funcao(x_atual + h0*dir_x, y_atual + h0*dir_y)
        	g1 = funcao(x_atual + h1*dir_x, y_atual + h1*dir_y)
        	g2 = funcao(x_atual + h2*dir_x, y_atual + h2*dir_y)
        	! Cálculo do numerador e denominador da fórmula do vértice
        	num = (h0**2 - h2**2)*g1 + (h2**2 - h1**2)*g0 + (h1**2 - h0**2)*g2
        	den = (h0 - h2)*g1 + (h2 - h1)*g0 + (h1 - h0)*g2
        	! Fallback
        	if (abs(den).lt.1.0e-12) then	! Caso o denominador seja muito próximo de zero,
        		passo = h0		! é feita a escolha do melhor passo entre h0, h1 e h2
        		g_max = g0
        		if (g1.gt.g_max) then
        			passo = h1
        			g_max = g1
        		end if
        		if (g2.gt.g_max) then
        			passo = h2
        			g_max = g2
        		end if
        	else
        		passo = 0.5d0*(num/den)
        	
        		if (passo.lt.0.0d0) then ! Garantia do passo não ser negativo
        			passo = h1
        		end if
        	end if
	end subroutine interpolacao
	
end program ppc4
