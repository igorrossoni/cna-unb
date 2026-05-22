# Configura o Gnuplot para modo de mapeamento de cores 2D
set pm3d map
set pm3d corners2color mean

# Escolha uma paleta de cores bonita (estilo 'jet' do azul ao vermelho)
set palette rgb 21,22,23

# Barra de cores
set colorbox
set cblabel "Temperatura [°C]"

# Configuração dos eixos
set title "Campo de Cores - Transferência de Calor na Parede do Reator com geração interna"
set xlabel "Posição x (m)"
set ylabel "Tempo (s)"
set grid

# Salva direto em uma imagem PNG
set terminal pngcairo size 800,600
set output 'campo_com_geracao.png'

# Plota usando a coluna 2 (x_pos) no eixo X e coluna 3 (Temperatura) como cor
splot "saida.dat" using 3:1:4 title ""
