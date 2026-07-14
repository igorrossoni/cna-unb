# PPC 6 - Transferência de Calor Bidimensional em Aleta
O programa `ppc6.f90` resolve numericamente a Equação de Laplace bidimensional aplicada a um problema de condução de calor em uma aleta retangular utilizando o Método das Diferenças Finitas. 

O programa funciona da seguinte maneira: 
* O usuário dá os parâmetros geométricos e físicos (comprimento, espessura, condutividade, coeficiente convectivo, temperaturas), as definições da malha (nós em x e y) e os parâmetros numéricos (tolerância de convergência e fator de relaxação);
* O programa discretiza o domínio gerando uma malha estruturada 2D e avalia as condições de contorno mistas: temperatura prescrita na base (Dirichlet) e troca convectiva de calor nas demais superfícies (Robin);
* O sistema linear resultante é resolvido utilizando três abordagens diferentes para comparação de desempenho: Eliminação de Gauss, Método de Liebmann (Gauss-Seidel) e Método de Liebmann com relaxação (SOR). Os métodos iterativos são otimizados operando diretamente sobre a malha de cálculo, dispensando o uso de matrizes densas. O processo iterativo se repete até que o erro seja menor que a tolerância definida;
* Após encontrar o campo de temperaturas convergido, os resultados nodais de cada método são exportados para arquivos `.dat`;
* O programa também isola a linha central da aleta ($y = H/2$) e calcula a solução analítica 1D clássica. Em seguida, extrai-se o erro percentual local relacionando a solução teórica 1D com a computacional 2D;
* Após a geração dos dados numéricos, utiliza-se o script `plot.gnu` para gerar as representações visuais através de mapas de calor e gráficos de linha 2D.

## Variáveis
**Parâmetros de entrada (dado pelo usuário)**
* `L` e `H`: comprimento e espessura da aleta, respectivamente;
* `cond_k`: condutividade térmica do material da aleta;
* `h_conv`: coeficiente convectivo de troca de calor com o ambiente;
* `Tb`: temperatura prescrita na base da aleta (parede fixa);
* `T_inf`: temperatura do fluido no ambiente livre;
* `nx` e `ny`: número de nós nas direções $x$ e $y$ da grade computacional;
* `tol`: tolerância máxima para convergência do erro nos métodos iterativos;
* `omega`: fator de sobre-relaxação ($\omega$) aplicado no método SOR.

**Controle iterativo e processamento**
* `A` e `b`: matriz de coeficientes e vetor de termos independentes, preenchidos exclusivamente para o método de Eliminação de Gauss;
* `iter_lieb` e `iter_relax`: contadores do número de iterações realizadas até atingir a convergência;
* `erro_lieb` e `erro_relax`: erro iterativo máximo (diferença entre o passo numérico atual e o anterior) obtido na matriz no ciclo atual;
* `tempo_gauss`, `tempo_lieb` e `tempo_relax`: tempo computacional de execução de cada subrotina de resolução matemática;
* `T_gauss`, `T_lieb` e `T_relax`: vetores/matrizes que armazenam o campo termal (solução do problema) obtido por cada um dos três métodos numéricos.

**Parâmetros físicos e adimensionais calculados**
* `n_total`: quantidade total de nós do domínio ($nx \times ny$);
* `dx` e `dy`: passo da malha (distância entre nós adjacentes) nas direções $x$ e $y$;
* `alpha`, `beta`, `gamma` e `phi`: constantes aglutinadas oriundas da discretização das condições de contorno de Robin (convecção) aplicadas nas bordas;
* `m_aleta`: constante do coeficiente global de transferência de calor da aleta aproximada para uma formulação bidimensional com profundidade unitária;
* `T_analitico`: temperatura prevista pela formulação teórica analítica 1D exata para o escoamento de calor em aletas retangulares;
* `erro_perc`: erro percentual comparativo entre o dado numérico (Bidimensional) e o dado analítico (Unidimensional).

## Inputs
* Parâmetros geométricos, físicos e numéricos digitados no terminal pelo usuário;
* Arquivos `gauss.dat`, `liebmann.dat`, `relax.dat` e `centro.dat` atuarão como inputs para a geração dos gráficos utilizando o Gnuplot.

## Outputs
* Os tempos de execução (`tempo`), número de iterações (`iter`) e erro de convergência (`erro`) de cada método computacional avaliado são impressos em tela no terminal de forma tabular;
* `gauss.dat`, `liebmann.dat`, `relax.dat`: arquivos numéricos contendo o campo completo. Organizados nas colunas: coordenada $x$, coordenada $y$ e Temperatura calculada $T(x,y)$;
* `centro.dat`: arquivo numérico da extração da linha média. Contém as colunas: $x$, $T_{numérico}$, $T_{analítico}$ e Erro Relativo (%), respectivamente;
* `mapa_temp.png`: imagem gerada pelo Gnuplot apresentando um painel comparativo dos 3 métodos via mapas de calor (superfícies de temperatura) e suas respectivas curvas de contorno (isotermas);
* `linha_central.png`: gráfico clássico 2D em linhas, justapondo e comparando o declínio de temperatura ao longo de $x$ entre o modelo numérico e o teórico;
* `erro_percentual.png`: gráfico ilustrando o comportamento espacial do erro percentual.

## Execução
A execução pode ser feita pelo terminal:
1. `gfortran ppc6.f90 -o ppc6`
2. `./ppc6`

    O usuário digita os parâmetros de entrada pedidos. Após o processamento numérico, os arquivos `.dat` são gerados na mesma pasta.

3. `gnuplot plot.gnu`

    Aqui são gerados os arquivos de imagem visuais (formato `.png`).

## Validação metodológica

### Parâmetros de Validação
* **Comprimento da aleta (`L`):** 1.0
* **Espessura da aleta (`H`):** 1.0
* **Condutividade térmica (`cond_k`):** 40.0
* **Coeficiente convectivo (`h_conv`):** 10.0
* **Temperatura da base (`Tb`):** 150.0
* **Temperatura ambiente (`T_inf`):** 25.0
* **Tolerância de convergência (`tol`):** 1.0d-5

### 1. Validação de Desempenho Computacional e Convergência

Para validar o desempenho dos algoritmos implementados, realizou-se uma análise comparativa do tempo de execução e do esforço computacional (número de iterações) exigido por cada método. Os testes foram conduzidos utilizando um domínio quadrado com L = 1 e H = 1, discretizado por uma malha espacial de 21x21 nós, sob uma tolerância de convergência rigorosa. Para o método de Liebmann com relaxação, fixou-se inicialmente o fator de sobre-relaxação em 1.5.

**Tabela 1:** Comparativo de desempenho entre os métodos de resolução (Malha 21x21; omega = 1.5).

| Método Matemático | Tempo de Execução (s) | Iterações | Erro Iterativo Final |
| :--- | :--- | :--- | :--- |
| Eliminação de Gauss | 0.01096 | - | - |
| Liebmann (Gauss-Seidel) | 0.02054 | 2524 | 9.9672E-06 |
| Liebmann com Relaxação | 0.00916 | 885 | 9.9055E-06 |

Os dados da Tabela 1 revelam uma dinâmica computacional muito interessante para esta resolução de malha (441 nós totais). O método iterativo clássico de Liebmann (Gauss-Seidel puro) apresentou o pior desempenho, sendo quase duas vezes mais lento que o método direto de Eliminação de Gauss (0.02054 s contra 0.01096 s). Isso ocorre porque, para malhas relativamente pequenas, o custo computacional de realizar milhares de iterações (2524 ciclos) acaba superando o custo das operações algébricas diretas de escalonamento e substituição da matriz.

Contudo, a introdução do fator de sobre-relaxação reverte esse cenário. A aplicação do método SOR com omega = 1.5 reduziu o número de iterações em cerca de 65% (de 2524 para 885), derrubando o tempo de execução para 0.00916 s e tornando a abordagem iterativa ligeiramente mais rápida que a eliminação direta.

Para investigar a fundo o ganho de eficiência atrelado ao método SOR, isolou-se o fator de relaxação (omega) como variável independente. Mantendo a mesma malha geométrica, o parâmetro foi variado de 1.0 (equivalente ao método de Liebmann clássico) até 1.9, aproximando-se do limite teórico de estabilidade para formulações de diferenças finitas.

**Tabela 2:** Efeito do fator de relaxação no esforço computacional e tempo de execução.

| Fator (omega) | Tempo de Execução (s) | Número de Iterações |
| :---: | :---: | :---: |
| 1.0 | 0.02078 | 2524 |
| 1.1 | 0.01841 | 2083 |
| 1.2 | 0.01520 | 1713 |
| 1.3 | 0.01360 | 1398 |
| 1.4 | 0.01038 | 1125 |
| 1.5 | 0.00697 |  885 |
| 1.6 | 0.00494 |  671 |
| 1.7 | 0.00405 |  477 |
| 1.8 | 0.00234 |  294 |
| 1.9 | 0.00141 |  160 |

A evolução dos dados comprova o enorme impacto da sobre-relaxação na aceleração da convergência, evidenciando um decaimento não linear estrito no número de iterações conforme o fator se aproxima de 2.0. Ao adotar omega = 1.9, o sistema convergiu em apenas 160 iterações, um esforço quase 16 vezes menor do que o exigido pelo método clássico.

Consequentemente, o tempo de execução sofreu uma redução dramática, caindo de 0.02078 s para apenas 0.00141 s. Essa análise de sensibilidade prova que, embora o método direto de Gauss seja bastante competitivo contra iterativos puros em malhas modestas, o método SOR devidamente calibrado é indiscutivelmente a ferramenta mais robusta e eficiente para a resolução computacional de malhas bidimensionais.

### 2. Estudo de Independência de Malha (Refinamento)

Para garantir que a solução numérica obtida não seja influenciada pelo erro de truncamento espacial da discretização, conduziu-se um estudo de refinamento de malha. O objetivo é demonstrar a convergência assintótica do campo de temperaturas e avaliar o impacto do tamanho dos elementos na precisão do modelo computacional.

Considerando os mesmos parâmetros de validação citados no início, avaliaram-se cinco níveis de densidade de malha, mantendo-se o fator de relaxação constante em $\omega = 1.5$. A métrica escolhida para atestar a estabilidade foi a temperatura no nó central da extremidade livre ($x = 1.0$, $y = 0.5$), acompanhada do erro percentual relativo à solução teórica analítica unidimensional (1D) nesse mesmo ponto.

**Tabela 3:** Estudo de convergência de malha avaliando a temperatura na extremidade ($x = 1.0$).

| Malha | Nós ($nx \times ny$) | $\Delta x = \Delta y$ | Temperatura em $x=L$ (°C) | Solução Analítica (°C) | Erro vs. Solução 1D (%) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | 11 $\times$ 11 | 0.1000 | 108.90 | 106.60 | 2.16 |
| 2 | 21 $\times$ 21 | 0.0500 | 108.87 | 106.60 | 2.13 |
| 3 | 41 $\times$ 41 | 0.0250 | 108.86 | 106.60 | 2.12 |
| 4 | 81 $\times$ 81 | 0.0125 | 108.84 | 106.60 | 2.11 |
| 5 | 101 $\times$ 101 | 0.0100 | 108.83 | 106.60 | 2.10 |

A observação da Tabela 3 comprova a convergência assintótica do modelo computacional. À medida que a malha é refinada, a temperatura calculada na extremidade estabiliza de forma suave, atingindo o patamar de 108.83 °C na malha mais densa. A variação da temperatura entre as malhas $41 \times 41$ e $101 \times 101$ é da ordem de centésimos de grau, comprovando que o modelo tornou-se independente da grade a partir de $\Delta x \le 0.025$. Por apresentar um excelente balanço entre baixo custo computacional e alta precisão assintótica, a malha de $41 \times 41$ nós é a configuração recomendada para esta geometria.

Além da validação numérica, os dados revelam um comportamento físico fundamental. O erro em relação à solução analítica 1D estabiliza na casa de 2.10%. Esse erro é substancialmente maior do que o observado em aletas finas tradicionais, e isso se deve puramente à proporção geométrica. Ao assumir $L = 1$ e $H = 1$ buscando facilitar os cálculos, a "aleta" deixa de ser delgada e passa a ser um bloco quadrado. 

Nessa condição de alta espessura, a premissa fundamental do modelo analítico 1D, de que a temperatura é uniforme ao longo de toda a seção transversal $y$, falha de maneira significativa. O intenso fluxo de calor convectivo nas bordas superior e inferior cria gradientes transversais severos que o modelo 1D ignora. A solução numérica 2D captura corretamente essa "resistência térmica" transversal, em que o núcleo da aleta se mantém mais quente do que as bordas, justificando o desvio de 2.10% como um retrato fiel da física bidimensional real do problema, e não como uma imprecisão computacional.


### 3. Validação Físico-Matemática e Discussão

Para a análise comparativa do campo térmico obtido pela simulação 2D com a solução analítica, utilizou-se o perfil de temperatura extraído ao longo da linha média transversal da aleta ($y = H/2$) proveniente do arquivo `centro.dat`, operando com a malha otimizada de $41 \times 41$ nós.

Os dados demonstram que ambas as soluções capturam o decaimento característico do problema, reduzindo a temperatura a partir da base em direção à extremidade livre sujeita à convecção. Contudo, o modelo numérico 2D prevê temperaturas consistentemente mais elevadas ao longo de toda a linha central em comparação ao modelo 1D. Na extremidade livre ($x = 1.0$), a formulação teórica prevê 106.60 °C, enquanto a simulação numérica crava 108.86 °C, estabelecendo um desvio relativo estabilizado em 2.12%. É possível verificar toda a distribuição de temperatura na linha central da aleta no gráfico da imagem `linha_central.png` e a distribuição de erros na `erro_percentual.png`.

Esse desvio percentual não representa um erro de discretização do algoritmo, mas sim o reflexo de uma limitação física da formulação analítica quando aplicada a geometrias espessas. O modelo analítico assume que a temperatura varia exclusivamente na direção longitudinal ($x$). Fisicamente, isso equivale a assumir que a aleta é infinitamente delgada transversalmente, de modo que a perda de calor por convecção seja sentida instantaneamente de forma uniforme por toda a seção transversal (ou seja, não há gradiente em $y$). A simulação 2D, regida pela Equação de Laplace acoplada às condições de Robin, não assume essa simplificação. Ela resolve o caminho percorrido pelo calor ao longo da espessura. Para que a energia escape do centro para o ambiente, ela precisa obrigatoriamente ser conduzida do núcleo ($y = 0.5$) até as superfícies limitantes ($y = 0$ e $y = 1.0$).

As condições de contorno de convecção impõem taxas reais de fluxo de calor nas bordas. Esse fenômeno gera gradientes térmicos verticais significativos: as superfícies da aleta, estando em contato direto com o fluido ambiente, tornam-se consideravelmente mais frias do que o núcleo do componente. 

No cenário atual de validação, a geometria consiste em um domínio de seção quadrada ($L = 1$ e $H = 1$). Por se tratar de um corpo muito espesso, a "resistência térmica condutiva" entre o centro e as bordas é elevada. A periferia da aleta acaba atuando como uma leve barreira isolante para o núcleo. Como o modelo 1D ignora essa resistência transversal, ele superestima a taxa de resfriamento na linha central, prevendo um decaimento térmico mais agressivo do que a realidade. Portanto, o desvio atestado de 2.12% quantifica a magnitude do efeito bidimensional de resistência térmica, provando que a simulação numérica possui uma fidelidade física superior para a predição termodinâmica de aletas espessas.

## Bibliografia
1. S. C. Chapra, R. P. Canale. “Métodos Numéricos para Engenharia.”, McGraw Hill, 5a edição (2008): 1-825.
2. R. G. Gontijo, “Notas de aula do curso de Cálculo Numérico Aplicado”, (2026).
3. Canal do YouTube Ciência e Brisa.