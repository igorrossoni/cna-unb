! PPC 6 - Transferência de Calor Bidimensional em Aleta
! Aluno: Igor de Oliveira Rossoni           Matrícula: 222031279
! A análise dos resultos está disponível na seção "Validação metodológica" do arquivo README.md

program ppc6
    implicit none
    
! Criação de variáveis
    integer, parameter :: iter_max = 100000
    integer :: nx, ny, n_total, i, j, k
    integer :: iter_lieb, iter_relax, j_centro
    real(8) :: L, H, cond_k, h_conv, Tb, T_inf
    real(8) :: dx, dy, beta, tol, omega
    real(8) :: alpha, gamma, phi
    real(8), allocatable :: A(:,:), b(:)
    real(8), allocatable :: T_gauss(:), T_lieb(:,:), T_relax(:,:), T_chute(:,:)
    real(8) :: tempo_gauss, tempo_lieb, tempo_relax
    real(8) :: erro_lieb, erro_relax
    real(8) :: x, y, m_aleta, theta, T_analitico, erro_perc, num, den

! Entrada de dados
    print *, 'Transferência de Calor em uma aleta'
    print *, 'Digite o comprimento da aleta L: '
    read *, L
    print *, 'Digite a espessura da aleta H: '
    read *, H
    print *, 'Digite a condutividade térmica k: '
    read *, cond_k
    print *, 'Digite o coeficiente convectivo h: '
    read *, h_conv
    print *, 'Digite a temperatura da base (parede esquerda) Tb: '
    read *, Tb
    print *, 'Digite a temperatura ambiente T_inf: '
    read *, T_inf
    print *, 'Digite o número de nós na direção x: '
    read *, nx
    print *, 'Digite o número de nós na direção y: '
    read *, ny
    print *, 'Digite a tolerância de convergência: '
    read *, tol
    print *, 'Digite o fator de relaxação omega: '
    read *, omega

! Parâmetros utilizados para a validação do problema (tirar do comentário caso queira replicar o resultado apresentado)
    ! L = 1.0d0
    ! H = 1.0d0
    ! cond_k = 40.0d0
    ! h_conv = 10.0d0
    ! Tb = 150.0d0
    ! T_inf = 25.0d0
    ! nx = 41
    ! ny = 41
    ! tol = 1.0d-5
    ! omega = 1.5d0

! Geração da malha
    n_total = nx*ny     ! Tamanho da matriz [A]
    dx = L/(nx - 1)     ! Passo na direção x
    dy = H/(ny - 1)     ! Passo na direção y

! Alocação das matrizes e vetores
    allocate(A(n_total,n_total))
    allocate(b(n_total))
    allocate(T_gauss(n_total))
    allocate(T_lieb(nx,ny))
    allocate(T_relax(nx,ny))
    allocate(T_chute(nx,ny))
    T_chute = T_inf

! Parâmetros da equação
    beta = h_conv*dx/cond_k     ! Número de Biot
    ! Os parâmetros alpha, gamma e phi servem para não sobrecarregar a geração da matriz [A]
    alpha = - 4.0d0
    gamma = - 2.0d0*(2.0d0 + beta)
    phi = - 4.0d0*(1.0d0 + beta)

! Montagem de [A] e {b}
    A = 0.0d0
    b = 0.0d0
    do j = 1, ny
        do i = 1, nx
            k = idx(i, j, nx)

            ! Base da aleta
            if (i.eq.1) then
                A(k,k) = 1.0d0
                b(k) = Tb

            ! Canto inferior direito
            else if (i.eq.nx.and.j.eq.1) then
                A(k, k) = phi
                A(k, idx(i-1, j, nx)) = 2.0d0 ! Nó da esquerda
                A(k, idx(i, j+1, nx)) = 2.0d0 ! Nó de cima
                b(k) = -4.0d0 * beta * T_inf

            ! Canto superior direito
            else if (i.eq.nx.and.j.eq.ny) then
                A(k, k) = phi
                A(k, idx(i-1, j, nx)) = 2.0d0  ! Nó da esquerda
                A(k, idx(i, j-1, nx)) = 2.0d0  ! Nó de baixo
                b(k) = -4.0d0 * beta * T_inf

            ! Borda inferior
            else if (j.eq.1) then
                A(k, k) = gamma
                A(k, idx(i-1, j, nx)) = 1.0d0  ! Nó da esquerda
                A(k, idx(i+1, j, nx)) = 1.0d0  ! Nó da direita
                A(k, idx(i, j+1, nx)) = 2.0d0  ! Nó de cima
                b(k) = -2.0d0 * beta * T_inf

            ! Borda superior
            else if (j.eq.ny) then
                A(k, k) = gamma
                A(k, idx(i-1, j, nx)) = 1.0d0  ! Nó da esquerda
                A(k, idx(i+1, j, nx)) = 1.0d0  ! Nó da direita
                A(k, idx(i, j-1, nx)) = 2.0d0  ! Nó de baixo
                b(k) = -2.0d0 * beta * T_inf

            ! Borda direita
            else if (i.eq.nx) then
                A(k, k) = gamma
                A(k, idx(i-1, j, nx)) = 2.0d0  ! Nó da esquerda
                A(k, idx(i, j-1, nx)) = 1.0d0  ! Nó de baixo
                A(k, idx(i, j+1, nx)) = 1.0d0  ! Nó de cima
                b(k) = -2.0d0 * beta * T_inf
            
            ! Nós internos
            else
                A(k, k) = alpha
                A(k, idx(i-1, j, nx)) = 1.0d0  ! Nó da esquerda
                A(k, idx(i+1, j, nx)) = 1.0d0  ! Nó da direita
                A(k, idx(i, j-1, nx)) = 1.0d0  ! Nó de baixo
                A(k, idx(i, j+1, nx)) = 1.0d0  ! Nó de cima
                b(k) = 0.0d0

            end if
        end do
    end do

! Resolução da matriz pelos 3 métodos
    write(*,'(A)') ""
    write(*,'(A)') "Realizando cálculos..."
    write(*,'(A)') "================================================================================="

    ! Método 1: Eliminação de Gauss
    call elim_gauss(n_total, A, b, T_gauss, tempo_gauss)
    write(*, '(A, F9.5, A)') "1. Gauss                  | Tempo: ", tempo_gauss, " s"
    
    ! Método 2: Liebmann (Gauss-Seidel)
    T_lieb = T_chute
    call liebmann(T_lieb, tol, iter_max, iter_lieb, erro_lieb, tempo_lieb)
    write(*, '(A, F9.5, A, I6, A, ES11.4)') "2. Liebmann               | Tempo: ", &
         tempo_lieb, " s | Iter: ", iter_lieb, " | Erro: ", erro_lieb
    
    ! Método 3: Liebmann com Relaxação
    T_relax = T_chute
    call lieb_relax(T_relax, tol, iter_max, iter_relax, erro_relax, tempo_relax)
    write(*, '(A, F9.5, A, I6, A, ES11.4)') "3. Liebmann com relaxação | Tempo: ", &
         tempo_relax, " s | Iter: ", iter_relax, " | Erro: ", erro_relax

    write(*,'(A)') "================================================================================="
    write(*,'(A)') ""

! Exportação dos dados para geração gráfica
    open(unit=10, file='gauss.dat', status='replace')
    open(unit=20, file='liebmann.dat', status='replace')
    open(unit=30, file='relax.dat', status='replace')
    do j = 1, ny
        do i = 1, nx
            k = idx(i, j, nx)
            x = (i-1)*dx
            y = (j-1)*dy
            write(10, *) x, y, T_gauss(k)
            write(20, *) x, y, T_lieb(i,j)
            write(30, *) x, y, T_relax(i,j)
        end do
        write(10, *) ""
        write(20, *) ""
        write(30, *) ""
    end do
    close(10)
    close(20)
    close(30)

! Distribuição de temperatura ao longo da linha média da aleta
    open(unit=40, file='centro.dat', status='replace')

    ! Para uma aleta 2D com espessura unitaria, o perimetro P = 2 e a Area Ac = H
    ! Logo m = sqrt(hP / kAc) = sqrt(2h / kH)
    m_aleta = sqrt((2.0d0*h_conv)/(cond_k*H))

    j_centro = (ny + 1)/2   ! Escolhe a linha mais proxima do centro (y = H/2)

    do i = 1, nx
        x = (i - 1)*dx
        ! Cálculo da solução analítica:
        num = cosh(m_aleta * (L - x)) + (h_conv / (m_aleta * cond_k)) * sinh(m_aleta * (L - x))
        den = cosh(m_aleta * L) + (h_conv / (m_aleta * cond_k)) * sinh(m_aleta * L)
        theta = num / den
        T_analitico = theta*(Tb - T_inf) + T_inf

        ! Cálculo do erro percentual
        ! É usado os valores de temperatura do método de Liebmann com sobre-relaxação, mas poderia ser utilizado qualquer um dos métodos
        erro_perc = abs((T_relax(i,j_centro) - T_analitico)/T_analitico)*100.0d0
        write(40, *) x, T_relax(i,j_centro), T_analitico, erro_perc
    end do
    close(40)

! Definição de funções e subrotinas
    contains
        ! Função inteira para indexação
        integer function idx(i, j, nx)
            integer, intent(in) :: i, j, nx
            idx = (j - 1) * nx + i
        end function idx

        ! Subrotina para aplicação do método de Eliminação Gaussiana
        subroutine elim_gauss(N, A_in, b_in, T, tempo)
            implicit none

            ! Variáveis de entrada e saída
            integer, intent(in) :: N                  ! Quantidade de nós
            real(8), intent(in) :: A_in(N,N), b_in(N) ! Matriz [A] e vetor de resultados {b}
            real(8), intent(out) :: T(N), tempo       ! Vetor {T} e tempo computacional total

            ! Variáveis locais da rotina
            real(8) :: A(N,N), b(N)
            real(8) :: fator, soma, t_inicio, t_final
            integer :: i, j, k

            call cpu_time(t_inicio)
            ! Substituição para não mudar a matriz original
            A = A_in
            b = b_in

            ! Eliminação progressiva
            do k = 1,n-1                    ! linha anterior
                do i=k+1,n                  ! linha atual
                    if (A(i,k).ne.0.0d0) then
                        fator=A(i,k)/A(k,k)     ! fator
                        do j=k+1,n              ! colunas
                            A(i,j) = A(i,j) - fator*A(k,j)  ! eliminando [A]
                        end do
                        b(i) = b(i) - fator*B(k)    ! alterando {b}
                    end if
                end do
            end do

            ! Substituição regressiva
            T(N) = b(N)/A(N,N)      ! començando de baixo
            do i=N-1,1, -1
                soma = b(i)         ! subindo de baixo para cima				
                do j= i+1,N
                    soma = soma - A(i,j)*T(j)   ! soma
                end do
                T(i) = soma/A(i,i)  ! substituindo em {T}
            end do

            call cpu_time(t_final)
            tempo = t_final - t_inicio

        end subroutine elim_gauss

        ! Subrotina para aplicação do método de Liebmann (Gauss-Seidel)
        subroutine liebmann(T, tol, iter_max, iter, erro_final, tempo)
            implicit none

            ! Variávies de entrada e saída
            real(8), intent(in) :: tol
            integer, intent(in) :: iter_max
            real(8), intent(inout) :: T(nx,ny)
            integer, intent(out) :: iter
            real(8), intent(out) :: erro_final, tempo

            ! Variáveis locais
            integer :: i, j
            real(8) :: soma, T_calc, erro
            real(8) :: t_inicio, t_final

            call cpu_time(t_inicio)
            iter = 0
            erro_final = 1.0d0

            do while (iter.lt.iter_max.and.erro_final.gt.tol)
                iter = iter + 1
                erro_final = 0.0d0
                ! Cálculo dos valores de temperatura em cada nó
                do j = 1, ny
                    do i = 1, nx
                        ! Base da aleta
                        if (i.eq.1) then
                            T_calc = Tb

                        ! Canto inferior direito
                        else if (i.eq.nx.and.j.eq.1) then
                            soma = 2.0d0*T(i-1,j) + 2.0d0*T(i,j+1)
                            T_calc = (- 4*beta*T_inf - soma)/phi

                        ! Canto superior direito
                        else if (i.eq.nx.and.j.eq.ny) then 
                            soma = 2.0d0*T(i-1,j) + 2.0d0*T(i,j-1)
                            T_calc = (-4.0d0*beta*T_inf - soma)/phi
                        
                        ! Borda inferior
                        else if (j.eq.1) then
                            soma = T(i-1,j) + T(i+1,j) + 2.0d0*T(i,j+1)
                            T_calc = (-2.0d0*beta*T_inf - soma)/gamma

                        ! Borda superior
                        else if (j.eq.ny) then
                            soma = T(i-1,j) + T(i+1,j) + 2.0d0*T(i,j-1)
                            T_calc = (-2.0d0*beta*T_inf - soma)/gamma
                        
                        ! Borda direita
                        else if (i.eq.nx) then
                            soma = T(i,j-1) + 2.0d0*T(i-1,j) + T(i,j+1)
                            T_calc = (-2.0d0*beta*T_inf - soma)/gamma

                        ! Nós internos
                        else
                            soma = T(i-1,j) + T(i+1,j) + T(i,j-1) + T(i,j+1)
                            T_calc = -soma/alpha
                        end if

                        erro = abs(T_calc - T(i,j))
                        if(erro.gt.erro_final) then
                            erro_final = erro
                        end if

                        T(i,j) = T_calc
                    end do
                end do
            end do

            call cpu_time(t_final)
            tempo = t_final - t_inicio
            
        end subroutine liebmann

        ! Subrotina para aplicação do método de Liebmann com relaxamento
        subroutine lieb_relax(T, tol, iter_max, iter, erro_final, tempo)
            implicit none

            ! Variávies de entrada e saída
            real(8), intent(in) :: tol
            integer, intent(in) :: iter_max
            real(8), intent(inout) :: T(nx,ny)
            integer, intent(out) :: iter
            real(8), intent(out) :: erro_final, tempo

            ! Variáveis locais
            integer :: i, j
            real(8) :: soma, T_calc, erro
            real(8) :: t_inicio, t_final

            call cpu_time(t_inicio)
            iter = 0
            erro_final = 1.0d0

            ! Cálculo dos valores de temperatura
            do while (iter.lt.iter_max.and.erro_final.gt.tol)
                iter = iter + 1
                erro_final = 0.0d0
                do j = 1, ny
                    do i = 1, nx
                        ! Base da aleta
                        if (i.eq.1) then
                            T_calc = Tb

                        ! Canto inferior direito
                        else if (i.eq.nx.and.j.eq.1) then
                            soma = 2.0d0*T(i-1,j) + 2.0d0*T(i,j+1)
                            T_calc = (- 4*beta*T_inf - soma)/phi

                        ! Canto superior direito
                        else if (i.eq.nx.and.j.eq.ny) then 
                            soma = 2.0d0*T(i-1,j) + 2.0d0*T(i,j-1)
                            T_calc = (-4.0d0*beta*T_inf - soma)/phi
                        
                        ! Borda inferior
                        else if (j.eq.1) then
                            soma = T(i-1,j) + T(i+1,j) + 2.0d0*T(i,j+1)
                            T_calc = (-2.0d0*beta*T_inf - soma)/gamma

                        ! Borda superior
                        else if (j.eq.ny) then
                            soma = T(i-1,j) + T(i+1,j) + 2.0d0*T(i,j-1)
                            T_calc = (-2.0d0*beta*T_inf - soma)/gamma
                        
                        ! Borda direita
                        else if (i.eq.nx) then
                            soma = T(i,j-1) + 2.0d0*T(i-1,j) + T(i,j+1)
                            T_calc = (-2.0d0*beta*T_inf - soma)/gamma

                        ! Nós internos
                        else
                            soma = T(i-1,j) + T(i+1,j) + T(i,j-1) + T(i,j+1)
                            T_calc = -soma/alpha
                        end if

                        erro = abs(T_calc - T(i,j))
                        if(erro.gt.erro_final) then
                            erro_final = erro
                        end if
                        
                        T(i,j) = omega*T_calc + (1.0d0 - omega)*T(i,j)
                    end do
                end do
            end do

            call cpu_time(t_final)
            tempo = t_final - t_inicio

        end subroutine lieb_relax

end program ppc6