# PPC 3 - Método de Thomas para resolver a transferência de calor em um reator
O programa "ppc3-thomas.f90" resolve o sistema linear [A]{T} = {B} para o problema de temperaturas em um reator utilizando o método de Thomas, descrito no arquivo "PPC3 e APC3.pdf". Além disso, serão plotados gráficos de cores e perfis para o problema considerando se o reator tem ou não geração interna.

## Variáveis
* L = espessura
* alpha = difusividade do material
* k_term = condutividade do material
* h_conv = coef. de transf. de calor por convecção
* T_inf = temperatura do fluido externo
* rho = densidade do material
* Cp = calor específico a pressão constante
* q_dot = taxa de geração
* A_fonte = termo fonte agrupado
* n = número de nós
* x_pos = posição na parede com relação ao número de nós
* dx = espaço entre nós
* dt = passo de tempo
* t_final = tempo total de simulação
* Fo = número de Fourier
* Bi = número de Biot
* e = vetor diagonal inferior
* f = vetor diagonal principal
* g = vetor diagonal superior
* r = vetor de fontes
* e_orig, f_orig, g_orig, r_orig = vetores originais
* T = vetor de temperatura no instante atual
* T_novo = vetor de temperatura no próximo passo de tempo
* p = passo de tempo atual
* p_total = número total de passos de tempo
* i e k = contadores inteiros

## Inputs
Para o método de Thomas, deve ser especificado os vetores e, f, g e r. Para o problema do reator, esses vetores já estão especificados diretamente no código, portanto, apenas os parâmetros físicos devem ser fornecidos. Para isso, deve-se alterar o valor das variáveis diretamente no código.

Para os cálculos sem geração interna, deve-se alterar o valor de q_dot para 0.0. Para calcular com geração interna, o valor utilizado para q_dot é 2.0e7. 

## Outputs
o programa "ppc3-thomas.f90" gera um arquivo "saida.dat" com os valores de p*dt (tempo em segundos), i (nó), x_pos (posição do nó na parede do reator em metros) e T(i) (temperatura para cada nó em °C) em colunas para cada bloco de tempo.

Para gerar os gráficos:
* **Caso sem geração interna:** usar os programas "campo_sem_geracao.gnu" e "evolucao_sem_geracao.gnu" para plotar o campo de cores ("campo_sem_geracao.png") e o perfil de temperaturas ("evolucao_sem_geracao.png");
* **Caso com geração interna:** usar os programas "campo_com_geracao.gnu" e "evolucao_com_geracao.gnu" para plotar o campo de cores ("campo_com_geracao.png") e o perfil de temperaturas ("evolucao_com_geracao.png");

## Execução
A execução desse conjunto de programas pode ser feito diretamente do terminal:
1. gfortran ppc3-thomas.f90
2. ./a.out

Aqui gera o arquivo "saida.dat"

* Caso sem geração interna:

3. gnuplot campo_sem_geracao.gnu evolucao_sem_geracao.gnu.

Aqui gera os arquivos campo_sem_geracao.png e evolucao_sem_geracao.png.

* Caso com geração interna:

3. gnuplot campo_com_geracao.gnu evolucao_com_geracao.gnu.

Aqui gera os arquivos campo_com_geracao.png e evolucao_com_geracao.png.

## Validação metodológica
* **Condições de contorno:** observando os gráficos da evolução de temperatura em função do tempo, a curva para o nó 7 responde instantaneamente no início do regime transiente, ela sai da temperatura inicial uniforme e segue para o equilíbrio com o fluido externo, validando o balanço de energia implícito no contorno. No caso sem geração interna, o domínio estabiliza de forma assintótica próximo à temperatura do fluido externo. A eliminação do gradiente térmico conforme o tempo passa mostra que o contorno esquerdo funciona como isolante térmico.

* **Geração interna:** para o caso sem geração, a temperatura máxima atinge cerca de 280 °C na superfície e o interior fica pouco abaixo pela condução pura a partir do contorno. Para o caso com geração, a temperatura máxima fica próximo de 800 °C, gerando o perfil parabólico característico em reatores com geração interna.

* **Estabilidade numérica:** a suavidade das curvas temporais e a distribuição dos nós intermediários confirmam a estabilidade numérica, pois o esquema implícito eliminou oscilações numéricas ou divergência com o passo de tempo adotado. Além disso, o algortimo de Thomas possui uma boa precisão numérica, resolvendo o sistema tridiagonal sem acúmulo de erros de arredondamento.

Os gráficos poderiam ter ficados ainda mais suaves, principalmente os campos de cores, se fossem utilizados mais nós e aumentado o tempo de simulação.

## Bibliografia
1. Roteiro PPC3 e APC3 disponibilizado pelo professor.
2. S. C. Chapra, R. P. Canale. “Métodos Numéricos para Engenharia.”, McGrawHill, 5a edição (2008): 1-825.
