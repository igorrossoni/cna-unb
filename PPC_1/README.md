# PPC 1 - Implementação computacional do método Runge-Kutta de quarta ordem
O objetivo desta atividade é implementar computacionalmente o método de Runge-Kutta de quarta ordem para a solução do problema de sedimentação de uma esfera em regime de baixo Reynolds.
O programa "ppc1-rk4.f90" é capaz de fazer a análise do problema de valor inicial para diferentes parâmetros, apenas alterando os valores de variáveis, como o número de Reynolds ou de Stokes, tanto de forma numérica quanto de forma analítica.

## Variáveis
* h = passo de tempo
* t0 e tf = tempos inicial e final
* t e v = variáveis adimensionais para o tempo e velocidade
* k1, k2, k3 e k4 = variáveis do método RK4
* St e Res = números de Stokes e Reynolds
* n = número de iterações
* v_analitica = solução analítica para o caso de Re = 0
* i = contador
* erro = erro absoluto para validação do ponto

## Inputs
Para cada análise, deve-se alterar o valor da variável diretamente no código, por exemplo:
* Caso de Re = 0, altera-se o valor de Res para 0 e varia-se o valor do número de Stokes. Além disso pode-se variar o passo de tempo h para analisar como o refinamento temporal afeta a qualidade da solução numérica.
* Caso de Re != 0, pode-se fazer uma análise alterando o valor de Res para comparação da velocidade para cada caso.

## Outputs
O programa "ppc1-rk4.f90' gera um arquivo "saida.dat" contendo os dados de tempo, velocidade, velocidade analítica e erro, respectivamente, no formato de colunas.

O programa "plot1.gnu" serve para plotar um gráfico para as diferentes análises a serem feitas, ele lê o arquivo "saida.dat" e plota as colunas necessárias no gráfico.

## Execução
A execução pelo terminal deve ser feita da seguinte maneira:
1. gfortran ppc1-rk4.f90
2. ./a.out
3. gnuplot plot1.gnu

## Validação metodológica
  Para o caso de Re = 0, observa-se que a solução numérica obtida pelo método RK4 coincide praticamente com a solução analítica ao longo de todo o intervalo de tempo analisado. As diferenças entre os valores aparece a partir da sexta casa decimal, sendo desprezível do ponto de vista prático. Graficamente, as curvas numérica e analítica se sobrepõem completamente.
  
  Para Re = 0, o sistema tende ao regime de Stokes, em que a solução é linear e a velocidade terminal é igual a 1. Conforme o número de Reynols aumenta, surge o termo não linear associado ao arrasto quadrático, diminuindo a velocidade terminal e afastando do caso linear. Assim, quanto maior o valor de Re, maior é o desvio em relação ao regime assintótico Re = 0.
  
  Para Re != 0, não há solução numérica simples em função do tempo, mas é possível determinar analiticamente a velocidade terminal impondo regime permanente, usando: v_analitica = (-1.0 + sqrt(1.0 + (3.0/2.0)*Res)) / ((3.0/4.0)*Res). A solução numérica obtida pelo método apresenta comportamento físico consistente, para Re = 0, a solução tende a 1, enquanto para Re != 0, converge para a velocidade terminal, que é menor que 1. Isso mostra confiabilidade no método RK4 e mostrando que ele consegue representar corretamente tanto o comportamento transitório quanto o regime permanente do sistema.
  
  O refinamento temporal, com a redução do passo de tempo, melhora significativamente a qualidade da solução numérica. Fazendo teste para h = 1, encontra-se um erro mais elevado, pois o passo é relativamente grande, reduzindo a capacidade do método de capturar a variação de velocidade ao longo do tempo, desse modo, a solução pode apresentar desvios em relação ao comportamento esperado. Ao reduzir o passo para h = 0.5, é perceptível a diferença na precisão, diminuindo o erro e, graficamente, tornando a curva mais suave. Com h = 0.1, o refinamento temporal é suficiente para garantir alta precisão, fazendo com que a solução numérica se aproxime muito da solução analítica, reduzindo o erro consideravelmente.

## Bibliografia
1. Roteiro PPC1 e APC1 disponibilizado pelo professor.
2. S. C. Chapra, R. P. Canale. “Métodos Numéricos para Engenharia.”, McGrawHill, 5a edição (2008): 1-825.
