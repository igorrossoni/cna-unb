# PPC 2 - Método de Bairstow para raízes de polinômios
O objetivo desta atividade é implementar computacionalmente o método de Bairstow para determinar as raízes de um polinômio, incluindo raízes reais e complexas. Além disso, o programa será aplicado para analisar o sistema dinâmico estudado na APC 2.

Dois programas foram escritos para esta atividade:
* "ppc2-bairstow.f90": a aplicação do método de Bairstow para encontrar raízes de polinômios de qualquer grau até 100.
* "varredura.f90": rotina que realiza uma varredura para valores de r0 e s0 no plano bidimensional, registrando o comportamento de convergência do método.

## Variáveis
### ppc2-bairstow.f90
* n = grau do polinômio
* a(0:100) = vetor que armazena os coeficientes do polinômio fornecido
* b(0:100) = vetor contendo os coeficientes da primeira divisão sintética do polinômio
* c(0:100) = vetor que recebe os coeficientes da segunda divisão sintética, utilizados para gerar o sistema linear do fator quadrático
* r e s = chutes iniciais e valores atualizados para os coeficientes do fator quadrático x²-rx+s
* dr e ds = correções calculados com base na regra de Cramer, servem para atualizar e refinar os valores de r e s
* det = calcula o determinante do sistema linear de c
* raizes(100) = vetor que armazena definitivamente os valores encontrados das raízes do polinômio (variável complexa)
* delta = discriminante usado para resolver a equação quadrática e extrair o par de raízes derivadas de r e s (variável complexa)
* i e k = contadores inteiros
* iter = contador de iterações
* total_iter = quantidade total de iterações
* max_iter = limite máximo de iterações, se extrapolar esse valor, o método diverge
* tol = tolência de erro permitido
* erro_r e erro_s = calculam o erro relativo associado à correção de r e s em cada passo

### varredura.f90
* n_inicial = grau do polinômio fornecido
* n = grau do polinômio atual
* a_inicial(0:100) = armazena permanentemente os coeficientes do polinômio fornecido, serve para resetar o polinômio de teste a cada ponto de varredura
* a(0:100) = vetor que armazena os coeficientes do polinômio fornecido, polinômio de teste
* b(0:100) = vetor contendo os coeficientes da primeira divisão sintética do polinômio
* c(0:100) = vetor que recebe os coeficientes da segunda divisão sintética, utilizados para gerar o sistema linear do fator quadrático
* r e s = variáveis que recebem o valor de r0 e s0 mas vao sendo atualizadas a cada iteração
* dr e ds = correções calculados com base na regra de Cramer, servem para atualizar e refinar os valores de r e s
* det = calcula o determinante do sistema linear de c
* tol = tolência de erro permitido
* erro_r e erro_s = calculam o erro relativo associado à correção de r e s em cada passo
* res = resolução da malha (500x500)
* i_grid e j_grid = contadores dos laços de varredura da malha, indo de 0 a 500
* r_min, r_max, s_min e s_max = definem os limites espaciais da varredura
* r0 e s0 = chute inicial na iteração atual da malha, são calculados proporcionalmente a i_grid e j_grid e nos limites máximo e mínimo
* i = contador inteiro
* iter = contador de iterações

## Inputs
* **ppc2-bairstow.f90**:

Ao rodar o programa, será necessário dar o input do grau n e de cada coeficiente a(i) do polinômio. Caso queira trocar os valores de chute inicial r e s, deve-se alterar no próprio código.

* **varredura.f90**:

Como este programa é feito para analisar o polinômio característico do sistema dinâmico desenvolvido na APC, não há inputs necessários. Caso queira fazer a varredura de algum outro polinômio, deve-se alterar os valores dos coeficientes a(i) diretamente no código.

## Outputs
* **ppc2-bairstow.f90**:

Saída direto do terminal mostrando a quantidade de iterações do método até encontrar as raízes e as raízes do polinômio descrito.

* **varredura.f90**:

Gera um arquivo "dados_varredura.csv" contendo r0, s0, r_final e s_final no formato de colunas.

* **gerar_fractal.plt**:

Pega os valores de r e s do arquivo "dados_varredura.csv" e plota um gráfico de cores conforme a convergência de cada valor chamado "fractal_bairstow.png".

## Execução
A execução pode ser feita direto do terminal a partir do passo a passo:
1. gfortran ppc2-bairstow.f90
2. ./a.out
3. inputs de n e a(i)

O output ocorre no próprio terminal

1. gfortran varredura.f90
2. ./a.out

Aqui é gerado o arquivo "dados_varredura.csv"

4. gnuplot gerar_fractal.plt

Aqui é gerado o arquivo "fractal_bairstow.png"

## Validação metodológica
O código "ppc2-bairstow.f90" foi validado usando o polinômio de grau 7:

P(x) = (x-1)(x-2)(x-3)(x-4)(x-5)(x-6)(x-7)

cujas raízes são: 1, 2, 3, 4 , 5, 6 e 7. Os valores encontrados de raízes para esse polinômio usando este código chegaram extremamente perto do verdadeiro, divergindo apenas por motivos de arredondamento.

O método é altamente sensível aos valores de chutes iniciais para r e s. Chutes mais próximos às raízes reais do sistema garantem uma convergência mais rápida. Chutes muito grandes podem causar instabilidade numérica, levando a um aumento na quantidade de iterações ou até mesmo à divergência.

Ao aplicar esse código ao polinômio característico do problema do massa-mola-amortecedor,

P(x) = 2x⁴ + 9x³ + 56x² + 70x + 200, encontra-se 4 raízes complexas:

x = -0.416 ± 2.20i e x = -1.833 ± 4.07i

A parte real negativa indica que o termo exponencial da solução decai com o tempo, portanto o sistema é estável e subamortecido. As partes imaginárias revelam as duas frequências naturais amortecidas do sistema acoplado, um sistema oscila lentamente (2.20 rad/s) e um rapidamente (4.07 rad/s) até que o movimento cesse completamente.

## Bibliografia
1. Roteiro PPC2 e APC2 disponibilizado pelo professor.
2. S. C. Chapra, R. P. Canale. “Métodos Numéricos para Engenharia.”, McGrawHill, 5a edição (2008): 1-825.
