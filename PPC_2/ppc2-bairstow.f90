! Igor de Oliveira Rossoni - 222031279

! Esse código funciona para qualquer polinômio de até grau 100 a partir dos valores de a(i) dados pelo usuário.

! Esse código foi validado usando o polinômio de grau 7:
! P(x) = (x-1)(x-2)(x-3)(x-4)(x-5)(x-6)(x-7)
! cujas raízes são: 1, 2, 3, 4 , 5, 6 e 7. Os valores encontrados de raízes para esse polinômio usando este código chegaram extremamente perto do verdadeiro, divergindo apenas por motivos de arredondamento.

! O método é altamente sensível aos valores de chutes iniciais para r e s. Chutes mais próximos às raízes reais do sistema garantem uma convergência mais rápida. Chutes muito grandes podem causar instabilidade numérica, levando a um aumento na quantidade de iterações ou até mesmo à divergência.

! Ao aplicar esse código ao polinômio característico do problema do massa-mola-amortecedor, P(x) = 2x⁴ + 9x³ + 56x² + 70x + 200, encontra-se 4 raízes complexas:
! x = -0.416 ± 2.20i
! x = -1.833 ± 4.07i
! A parte real negativa indica que o termo exponencial da solução decai com o tempo, portanto o sistema é estável e subamortecido. As partes imaginárias revelam as duas frequências naturais amortecidas do sistema acoplado, um sistema oscila lentamente (2.20 rad/s) e um rapidamente (4.07 rad/s) até que o movimento cesse completamente.

program bairstow
implicit none

integer n, i, iter, k, max_iter, total_iter
real a(0:100), b(0:100), c(0:100)
complex raizes(100), delta
real r, s, dr, ds, det
real tol, erro_r, erro_s

tol = 1.0e-6		! Tolerância
max_iter = 100		! Quantidade máxima de iterações, se passar, o método diverge
total_iter = 0		! Total de iterações

! Definição do polinômio pelo usuário
write(*,*) 'Digite o grau do polinômio'
read(*,*) n
do i = n, 0, -1
	write(*,'(A,I0,A)') 'Digite o coeficiente a(',i,'):'
	read(*,*) a(i)
end do

k = 0			! Contador de raízes

do while (n > 0)
! Enquanto n - 2 >= 3, utiliza-se o método de Bairstow
	if (n.ge.3) then
		r = 0.5
		s = 0.5
		do iter = 1, max_iter
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
			
! Caso o determinante do sistema seja igual a 0 ou muito pequeno, são dados novos valores para r e s e recomeça o laço
			if (abs(det).lt.1.0e-7) then
				r = r + 0.5
				s = s + 0.5
				cycle
			end if
! Definição de dr e ds pelo método de Cramer dos determinantes
			dr = (-(-b(1)*c(2)) + (-b(0)*c(3)))/det
			ds = (-(-b(0)*c(2)) + (-b(1)*c(1)))/det

! Atualizando os valores de r e s
			r = r + dr
			s = s + ds
! Cálculo do erro relativo
			erro_r = abs(dr/r)
			erro_s = abs(ds/s)
			
! Se o erro relativo for menor que a tolerância, o laço termina, se for maior, ele recomeça
			if (erro_r.lt.tol.and.erro_s.lt.tol) exit
		end do
		total_iter = total_iter + iter
	
! Resolver a eq. quadrática
		delta = cmplx(r*r + 4.0*s, 0.0)
		raizes(k+1) = (cmplx(r,0.0) + sqrt(delta))/2.0
		raizes(k+2) = (cmplx(r,0.0) - sqrt(delta))/2.0
		k = k + 2

! Deflação polinomial
		do i = 0, n-2
			a(i) = b(i+2)
		end do
		n = n - 2

! Caso n - 2 = 2, resolve-se por Bhaskara
	else if (n == 2) then
		r = a(1)/a(2)
		s = a(0)/a(2)
		delta = cmplx(r*r - 4.0*s, 0.0)
		raizes(k+1) = (-cmplx(r,0.0) + sqrt(delta))/2.0
		raizes(k+2) = (-cmplx(r,0.0) - sqrt(delta))/2.0
		k = k + 2
		
		n = n - 2
	
! Caso n - 2 = 1, x_r = -s/r
	else if (n ==1) then
		raizes(k+1) = -a(0)/a(1)
		k = k + 1
		
		n = n - 1
		
	end if
end do

write(*,'(A, I0)') 'Total de iterações:', total_iter
write(*,*) 'As raízes do polinômio são:'
do i = 1, k
	write(*,*) raizes(i)
end do

end program bairstow

