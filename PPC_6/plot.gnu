# ==========================================================
# CONFIGURAÇÃO GERAL E INÍCIO DO MULTIPLOT
# ==========================================================
set terminal pngcairo size 1800,500 enhanced font 'Arial,12'
set output 'mapa_temp.png'

# Ativa o multiplot com 1 linha e 3 colunas
set multiplot layout 1,3 title "Análise Térmica da Aleta (malha 41x41)" font ",16"

# ==========================================================
# GRÁFICO 1: Eliminação Gaussiana
# ==========================================================
set title "Eliminação Gaussiana" font ",14"
set xlabel "Comprimento L"
set ylabel "Espessura H"

set view map
set contour base
set cntrparam levels auto 12
unset surface

# Paleta de cores e configuracoes 3D
set palette defined (0 "dark-blue", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red", 6 "dark-red")
set pm3d at b # interpolate 2,2

splot 'gauss.dat' using 1:2:3 with pm3d notitle, \
      'gauss.dat' using 1:2:3 with lines nosurface linecolor rgb "black" notitle

# ==========================================================
# GRÁFICO 2: Liebmann
# ==========================================================
set title "Método de Liebmann" font ",14"
set xlabel "Comprimento L"
set ylabel "Espessura H"

set view map
set contour base
set cntrparam levels auto 12
unset surface

# Paleta de cores e configuracoes 3D
set palette defined (0 "dark-blue", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red", 6 "dark-red")
set pm3d at b # interpolate 2,2

splot 'liebmann.dat' using 1:2:3 with pm3d notitle, \
      'liebmann.dat' using 1:2:3 with lines nosurface linecolor rgb "black" notitle

# ==========================================================
# GRÁFICO 3: Liebmann com relaxação
# ==========================================================
set title "Método de Liebmann com relaxação" font ",14"
set xlabel "Comprimento L"
set ylabel "Espessura H"

set view map
set contour base
set cntrparam levels auto 12
unset surface

# Paleta de cores e configuracoes 3D
set palette defined (0 "dark-blue", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red", 6 "dark-red")
set pm3d at b # interpolate 2,2

splot 'relax.dat' using 1:2:3 with pm3d notitle, \
      'relax.dat' using 1:2:3 with lines nosurface linecolor rgb "black" notitle

unset multiplot

# ==========================================================
# GRÁFICO 4: Comparação na Linha Central (Novo Arquivo)
# ==========================================================
# Redefine o tamanho para um gráfico simples e nomeia o novo arquivo
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'linha_central.png'

# Limpa as configurações de mapa de calor (3D) para voltar ao 2D
unset pm3d
unset contour
unset view
set surface

# Configuração de eixos, título e legenda
set title "Distribuição de Temperatura na Linha Central (y = H/2) (malha 41x41)" font ",14"
set xlabel "Posição x na aleta"
set ylabel "Temperatura"
set grid
set key top right box

# Plota os dados bidimensionais de centro.dat
# Coluna 1 = x, Coluna 2 = T numérico, Coluna 3 = T analítico
plot 'centro.dat' using 1:2 with linespoints pt 7 ps 1.2 linecolor rgb "blue" title "Numérico (Relaxação)", \
     'centro.dat' using 1:3 with lines lw 2 linecolor rgb "red" title "Solução Analítica"

# ==========================================================
# GRÁFICO 5: Comportamento do Erro Percentual
# ==========================================================
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'erro_percentual.png'

set title "Erro Relativo entre Modelo Numérico 2D e Analítico 1D (malha 41x41)" font ",14"
set xlabel "Posição x na aleta"
set ylabel "Erro (%)"
set grid

# Plota a coluna 1 (x) contra a coluna 4 (erro)
plot 'centro.dat' using 1:4 with linespoints pt 7 ps 1.2 lw 2 linecolor rgb "forest-green" title "Erro (%)"