# Configurações gerais de estilo e saída
set terminal pngcairo size 1200, 500 enhanced font 'Verdana,11'
set output 'perfis.png'

# Define que o gráfico será dividido em 2 subplots (1 linha, 2 colunas)
set multiplot layout 1, 3 title "Perfis de similaridade"

set grid
set key bottom right

# No arquivo, Coluna 1 = eta, Coluna 2 = y1, Coluna 3 = y2, Coluna 4 = y3
# -----------------------------------------------------------------
# GRÁFICO 1: f(η)
# -----------------------------------------------------------------
set xlabel "η"
set ylabel "y1"

plot 'saida.dat' using 1:2 with lines title "f(η)" ls 1

# -----------------------------------------------------------------
# GRÁFICO 2: f'(η)
# -----------------------------------------------------------------
set xlabel "η"
set ylabel "y2"

plot 'saida.dat' using 1:3 with lines title "f'(η)" ls 1

# -----------------------------------------------------------------
# GRÁFICO 3: f''(η)
# -----------------------------------------------------------------
set xlabel "η"
set ylabel "y3"

plot 'saida.dat' using 1:4 with lines title "f''(η)" ls 1

unset multiplot