# PPC 4 - Maximização de duas variáveis
O programa "ppc4.f90" resolve o problema de maximização de duas variáveis utilizando os métodos do Aclive Máximo e dos Gradientes Conjugados (Fletcher-Reeves). A função a ser maximizada é f(x,y) = 2xy + 2x - x² - 2y², cujo ótimo analítico é (x, y) = (2, 1). O código foi validado com os pontos iniciais (x0, y0) = (-2, 3), mas funciona para qualquer ponto inicial de busca.

O programa "grafico.py" é um código em python para a geração do gráfico das curvas de nível com a comparação entre os dois métodos, ele utiliza as bibliotecas numpy e matplotlib. Este código usa os três arquivos gerados pelo "ppc4.f90" para gerar o gráfico.

## Variáveis
* x, y = par ordenado
* x_inicial, y_inicial = par ordenado inicial de busca
* dx, dy = componentes do gradiente
* dx_ant, dy_ant = componentes do gradiente no passo anterior
* f_val = função avaliada no par ordenado
* passo = passo de distância
* erro = módulo do gradiente
* px, py = direções de busca (método dos gradientes conjugados)
* beta = parâmetro minimizador
* x_malha, y_malha, f_malha = parâmetros para a geração da malha
* h0, h1, h2, g0, g1, g2, num, den, g_max = parâmetros para interpolação quadrática

## Inputs
Ao rodar o ppc4.f90, o usuário deverá dar o input dos pontos iniciais de busca (x, y) pelo terminal. Para a validação do código, foi utilizado (x, y) = (-2, 3).

## Outputs
Durante a execução do ppc4.f90, é impresso na tela: iter, x, y, f_val, passo, erro. Além disso, gera três arquivos, "output1.dat" e "output2.dat" (logs dos métodos AM e GC, respectivamente) e "function.dat" (amostras de (x, y, f) para plotagem das curvas de nível).

O formato dos arquivos "output1.dat" e "output2.dat" é por colunas na seguinte formatação: iter, erro, passo, x, y, dx, dy. Já o "function.dat" é: x_malha, y_malha, f_malha.

O grafico.py gera um grafico "curvas_de_nivel.png", que mostra as curvas de nível da função e uma comparação entre as buscas dos dois métodos.

## Execução
A execução pode ser realizada diretamente pelo terminal:
1. gfortran ppc4.f90
2. ./a.out
   
   Aqui é gerado "output1.dat", "output2.dat" e "function.dat"
   
4. python grafico.py
   
   Aqui é gerado "curvas_de_nivel.png"

## Validação metodológica
O método do aclive máximo utiliza a direção pelo gradiente da função enquanto o método dos gradientes conjugados atualiza a direção de busca pelo gradiente da função no ponto atual somado a um parâmetro minimizador multiplicado pelo passo de busca anterior. O parâmetro minimizador é dado pela norma do gradiente atual ao quadrado divido pela norma do gradiente anterior ao quadrado.

Para calcular o passo dos métodos, foi utilizada interpolação quadrática na linha de busca, utilizando a fórmula do vértice, e, caso a parábola for degenerada, é adotado como passo o melhor entre h0, h1 e h2.

Fazendo a análise do gráfico, percebe-se que o método AC demora muito mais para convergir do que o método GC. Além disso, percebe-se a diferença na trajetória, com o AM sendo em zigue-zague e o GC sendo mais direto.

## Bibliografia
1. Roteiro PPC4 e APC4 disponibilizado pelo professor.
2. S. C. Chapra, R. P. Canale. “Métodos Numéricos para Engenharia.”, McGrawHill, 5a edição (2008): 1-825.
