! PPC 5 - Solução numérica para a Equação de Blasius usando o Método de Runge-Kutta de quarta ordem e o Método do Tiro
! Aluno: Igor de Oliveira Rossoni          Matrícula: 222031279

! A análise dos parâmetros encontrados está no arquivo README.md na parte de validação metodológica.
! GitHub: https://github.com/igorrossoni/cna-unb.git

program blasius
    implicit none
    integer :: i, n_passos
    logical :: gravar_arquivo = .false., achou_eta_99 = .false.

    ! Parâmetros dados pelo usuário
    real(8) :: s0, delta_eta, eta_max, tol
    integer :: max_iter

    ! Variáveis do Método do Tiro
    real(8) :: s_atual, s_new, s_old
    real(8) :: erro, erro_old
    integer :: iter

    ! Variáveis do Método RK4
    real(8), dimension(3) :: y, y_new
    real(8), dimension(3) :: k1, k2, k3, k4
    real(8) :: eta, eta_99, Cf, C_delta

    ! Leitura dos parâmetros iniciais
    write(*,*) '======================================================'
    write(*,*) '        CONFIGURAÇÃO DOS PARÂMETROS DE ENTRADA        '
    write(*,*) '======================================================'
    
    write(*, '(A)', advance='no') ' Chute inicial para f"(0) (ex: 1.0)       : '
    read *, s0
    
    write(*, '(A)', advance='no') ' Passo de integração delta_eta (ex: 0.01) : '
    read *, delta_eta
    
    write(*, '(A)', advance='no') ' Valor máximo de eta (ex: 20.0)           : '
    read *, eta_max
    
    write(*, '(A)', advance='no') ' Tolerância de convergência (ex: 1d-6)    : '
    read *, tol
    
    write(*, '(A)', advance='no') ' Número máximo de iterações (ex: 50)      : '
    read *, max_iter
    
    write(*,*) '======================================================'
    write(*,*) 'Aplicando o Método de Runge-Kutta de ordem 4 e o Método do Tiro.'
    write(*,*) ''

    n_passos = int(eta_max/delta_eta) + 1   ! Número de passos
    s_atual = s0                            ! Definindo o chute inicial como s0
    iter = 0
    erro = 100.0d0

    do while (abs(erro).gt.tol)

        iter = iter + 1
        if (iter.ge.max_iter) then
            print *, 'Não convergiu dentro do máximo de iterações!'
            exit
        end if

        call rk4
        erro = y(2) - 1.0d0 ! Cálculo do erro. Vale a pena preservar o sinal do erro para saber se está acima ou abaixo do valor de convergência

        ! Método da secante para atualização do chute
        if (iter.eq.1) then
            s_new = s_atual + 0.1d0
            s_old = s_atual
            erro_old = erro
        else
            s_new = s_atual - erro*(s_atual - s_old)/(erro - erro_old)
            s_old = s_atual
            erro_old = erro
        end if
        s_atual = s_new
        
    end do

    ! Gerar o arquivo de saída
    gravar_arquivo = .true.
    call rk4
    gravar_arquivo = .false.

    Cf = 2.0d0*s_atual      ! f"(0) = s_atual
    C_delta = eta_99        ! C_delta calculado numéricamente é igual ao eta_99

    ! Impressão dos resultados no terminal
    write(*,*) 'Arquivo saída.dat gerado.'
    write(*,*) '======================================================'
    write(*,*) '       RESULTADOS DA SIMULAÇÃO - CAMADA LIMITE        '
    write(*,*) '======================================================'
    
    write(*, '(A, F12.6)') ' Valor convergido de f"(0)          : ', s_atual
    write(*, '(A, I12)')   ' Número de iterações do Tiro        : ', iter
    write(*, '(A, E12.4)') ' Erro final na borda [f''(inf) - 1]  : ', erro
    write(*, '(A, F12.6)') ' Valor final alcançado f''(eta_max)  : ', y(2)
    
    write(*,*) '------------------------------------------------------'
    write(*,*) 'PARÂMETROS FÍSICOS (ADIMENSIONAIS)'
    write(*, '(A, F12.6)') ' Constante de Atrito [Cf*sqrt(Rex)] : ', Cf
    write(*, '(A, F12.2)') ' Coeficiente de Espessura [C_delta] : ', C_delta
    write(*,*) '======================================================'
    write(*,*) 'Para gerar os gráficos, utilize: gnuplot plot.gnu'

    contains
    subroutine rk4
        ! Atribuição dos valores iniciais, com f" igual ao chute inicial
        y(1) = 0.0d0
        y(2) = 0.0d0
        y(3) = s_atual
        eta = 0.0d0

        if (gravar_arquivo) then        ! Gravação do arquivo de resultados para o chute convergido, se gravar_arquivo = .true.
            open(unit=10, file='saida.dat', status='replace')
            write(10,*) eta, y(1), y(2), y(3)
        end if
        
        ! Runge-Kutta de 4° ordem
        do i=1, n_passos
            ! Cálculo de k1 para y1, y2 e y3
            k1(1) = y(2)
            k1(2) = y(3)
            k1(3) = - y(1)*y(3)/2.0d0

            ! Cálculo de k2 para y1, y2 e y3
            k2(1) = y(2) + delta_eta*(k1(2))/2.0d0
            k2(2) = y(3) + delta_eta*(k1(3))/2.0d0
            k2(3) = - (y(1) + delta_eta*k1(1)/2.0d0)*(y(3) + delta_eta*k1(3)/2.0d0)/2.0d0

            ! Cálculo de k3 para y1, y2 e y3
            k3(1) = y(2) + delta_eta*(k2(2))/2.0d0
            k3(2) = y(3) + delta_eta*(k2(3))/2.0d0
            k3(3) = - (y(1) + delta_eta*k2(1)/2.0d0)*(y(3) + delta_eta*k2(3)/2.0d0)/2.0d0

            ! Cálculo de k4 para y1, y2 e y3
            k4(1) = y(2) + delta_eta*k3(2)
            k4(2) = y(3) + delta_eta*k3(3)
            k4(3) = - (y(1) + delta_eta*k3(1))*(y(3) + delta_eta*k3(3))/2.0d0

            ! Cálculo do próximo valor de y1, y2 e y3
            y_new(1) = y(1) + (delta_eta/6.0d0)*(k1(1) + 2.0d0*k2(1) + 2.0d0*k3(1) + k4(1))
            y_new(2) = y(2) + (delta_eta/6.0d0)*(k1(2) + 2.0d0*k2(2) + 2.0d0*k3(2) + k4(2))
            y_new(3) = y(3) + (delta_eta/6.0d0)*(k1(3) + 2.0d0*k2(3) + 2.0d0*k3(3) + k4(3))
            y = y_new

            ! Próximo valor de eta
            eta = eta + delta_eta
            
            if (gravar_arquivo) then
                write(10,*) eta, y(1), y(2), y(3)
                if (y(2).ge.0.99d0.and..not.achou_eta_99) then  ! Usado para encontrar eta_99
                    eta_99 = eta
                    achou_eta_99 = .true.
                end if
            end if
        end do
        if (gravar_arquivo) close(10)
    end subroutine

end program blasius