! Igor de Oliveira Rossoni - 222031279

! Este programa faz a varredura de valores de r0 e s0 para a geração do fractal

program varredura
implicit none

! Parâmetros da geração da grade de análise
integer i_grid, j_grid
integer, parameter :: res = 500
real r0, s0
real, parameter :: r_min = -5.0, r_max = 5.0
real, parameter :: s_min = -5.0, s_max = 5.0

! Parâmetros do método de Bairstow
real a(0:100), b(0:100), c(0:100), a_inicial(0:100)
real r, s, dr, ds, det
real tol, erro_r, erro_s
integer n, n_inicial, i, iter

! Polinômio
n_inicial = 6
a_inicial(6) = -2
a_inicial(5) = 8
a_inicial(4) = -12
a_inicial(3) = -1067
a_inicial(2) = 23
a_inicial(1) = 4
a_inicial(0) = 53

tol = 1.0e-6

! Abertura do arquivo de saída de dados
open(10, file='dados_varredura.csv', status='replace')
write(10,*) "r0,s0,r_final,s_final"

do i_grid = 0, res
	r0 = r_min + (r_max - r_min) * real(i_grid)/real(res)
	do j_grid = 0, res
		s0 = s_min + (s_max - s_min) * real(j_grid)/real(res)
		
		! Reset de variáveis
		n = n_inicial
		a = a_inicial
		r = r0
		s = s0
		erro_r = 1.0
		erro_s = 1.0
		do iter = 0, 100
			! Cálculo dos coeficientes b
			b(n) = a(n)
			b(n-1) = a(n-1) + r*b(n)
			do i = n-2, 0, -1
				b(i) = a(i) + r*b(i+1) + s*b(i+2)
			end do
		
			! Cálculo dos coeficientes c
			c(n) = b(n)
			c(n-1) = b(n-1) + r*c(n)
			do i = n-2, 1, -1
				c(i) = b(i) + r*c(i+1) + s*c(i+2)
			end do
		
			! Sistema linear de c's e deltas:
			! c1*dr + c2*ds = -b0
			! c2*dr + c3*ds = -b1
			det = c(3)*c(1) - c(2)*c(2)
			
			! Definição de dr e ds pelo método de Cramer dos determinantes
			dr = (-(-b(1)*c(2)) + (-b(0)*c(3)))/det
			ds = (-(-b(0)*c(2)) + (-b(1)*c(1)))/det

			! Atualizando os valores de r e s
			r = r + dr
			s = s + ds
			! Cálculo do erro relativo
			erro_r = abs(dr/r)
			erro_s = abs(ds/s)
			if (abs(erro_r).lt.tol.and.abs(erro_s).lt.tol) exit
		end do
		
		! Escrita dos resultados encontrados
		write(10, '(4(F16.8, ","))') r0, s0, r, s
	end do
	
end do

close(10)
write(*,*) 'Arquivo dados_varredura.csv gerado'

end program varredura
