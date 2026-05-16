# Configurações de saída
set terminal pngcairo size 800,800 enhanced font 'Verdana,10'
set output 'fractal_bairstow.png'

# Configurações do gráfico
set title "Fractal de Bairstow"
set xlabel "r_0"
set ylabel "s_0"

# Define o separador como vírgula (para ler o CSV)
set datafile separator ","

# Paleta de cores (pode ajustar para 'set palette rgbformulae 33,13,10' para cores clássicas)
set palette rgbformulae 33,13,10

# O truque: usamos a soma ou combinação de rfim e sfim para definir a cor
# rfim está na coluna 3 e sfim na coluna 4
plot 'dados_varredura.csv' using 1:2:($3 + $4) with points pt 5 ps 1.0 lc palette notitle
