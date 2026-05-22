set terminal pngcairo size 1200,800
set output 'evolucao_com_geracao.png'

set title "Temperatura em função do tempo para cada nó - com geração"
set xlabel "Tempo (s)"
set ylabel "Temperatura (°C)"
set key bottom right

set grid

plot \
'saida.dat' using 1:($2==1 ? $4 : 1/0) with linespoints title '0.000', \
'saida.dat' using 1:($2==2 ? $4 : 1/0) with linespoints title '0.002', \
'saida.dat' using 1:($2==3 ? $4 : 1/0) with linespoints title '0.004', \
'saida.dat' using 1:($2==4 ? $4 : 1/0) with linespoints title '0.006', \
'saida.dat' using 1:($2==5 ? $4 : 1/0) with linespoints title '0.008', \
'saida.dat' using 1:($2==6 ? $4 : 1/0) with linespoints title '0.010', \
'saida.dat' using 1:($2==7 ? $4 : 1/0) with linespoints title '0.012'

set output
