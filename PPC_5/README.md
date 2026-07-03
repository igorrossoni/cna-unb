# PPC 5 - Equação de Blasius
O programa `blasius.f90` resolve numéricamente a Equação de Blasius utilizando o Método de Runge-Kutta de quarta ordem e o Método do Tiro.

O programa funciona da seguinte maneira: 
* O usuário dá os parâmetros de chute inicial, passo de integração, $\eta$ máximo, tolerância de convergência e número máximo de iterações do método do tiro;
* O programa aplica o RK4 para o valor do chute inicial até o valor de eta_max, calcula o erro e encontra um novo valor de chute usando o método da secante. O processo é repetido até o erro ser menor que a tolerância definida;
* Após encontrar o valor do chute que convergiu, é aplicado novamente o RK4 para gerar o arquivo "saida.dat" com os resultados de $\eta$, $f(\eta)$, $f'(\eta)$ e $f"(\eta)$. Durante essa aplicação do RK4, é encontrado $f'(\eta_{99})$ = 0.99, que corresponde ao coeficiente de espessura C_delta;
* Com o valor do chute convergido, encontra-se também o coeficiente local de atrito na parede (adimensional): $Cf = 2*f"(0)$.

Após a geração do arquivo "saida.dat", utiliza-se o programa `plot.gnu` para gerar um gráfico com os perfis de similaridade:
- $f(\eta)$;
- $f'(\eta)$;
- $f"(\eta)$.

## Variáveis
**Parâmetros de entrada (dado pelo usuário)**
* `s0`: chute inicial para a derivada segunda na parede;
* `delta_eta`: passo de integração;
* `eta_max`: limite do domínio de integração (condição longe da placa);
* `tol`: tolerância máxima para convergência do erro na condição de contorno;
* `max_iter`: quantidade máxima de iterações para o método do tiro.

**Controle iterativo e integração**
* `s_atual`, `s_new` e `s_old`: variáveis empregadas no método da secante para encontrar o próximo chute;
* `erro` e `erro_old`: diferença quantitativa calculada entre o limite teórico da velocidade no infito e o valor encontrado pela simulação no passe atual (dado por: $f'(eta_max) - 1$);
* `iter`: contador de iterações do método do tiro;
* `y` e `y_new`: vetores de dimensão 3, em que a primeira posição representa $f(\eta)$, a segunda $f'(\eta)$ e a terceira $f"(\eta)$;
* `k1`, `k2`, `k3` e `k4`: parâmetros usados no método de Runge-Kutta;
* `eta`: variável de similaridade adimensional.

**Parâmetros físicos adimensionais calculados**
* `eta_99`: posição adimensional em que a velocidade do perfil atinge 99% da velocidade do escoamento livre do fluido;
* `Cf`: constante do coeficiente de atrito local na parede;
* `C_delta`: coeficiente que dita o crescimento da espessura da camada limite.


## Inputs
* Parâmetros de entrada digitados no terminal;
* `saida.dat` para a geração dos gráficos usando o Gnuplot.

## Outputs
* `s_atual`, `iter`, `erro`, `y(2)`, `Cf` e `C_delta` são impressos no terminal;
* `saida.dat`: arquivo numérico gerado após a convergência do chute. Contém as colunas de dados: $\eta$, $y(1)$, $y(2)$ e $y(3)$, respectivamente;
* `perfis.png`: arquivo de imagem gerado pelo script do Gnuplot, exibindo visualmente as curvas em função de corrente e suas derivadas.

## Execução
A execução pode ser feita pelo terminal:
1. gfortran blasius.f90
2. ./a.out

    O usuário digita os parâmetros de entrada pedidos. Após o processamento numérico, é gerado o arquivo `saida.dat`.

3. gnuplot plot.gnu

    Aqui é gerado "perfis.png"

## Validação metodológica
Para a validação do programa e geração dos gráficos, foram utilizados os seguintes valores para os parâmetros de entrada:
- `s0` = 1.0
- `delta_eta` = 0.01
- `eta_max` = 20.0
- `tol` = 1d-6
- `max_iter` = 50

Os resultados encontrados foram:
- `f"(0)` = 0.332057
- `iter` = 7
- `erro` = 0.1161E-06
- `f'(eta_max)` = 1.000000
- `Cf` = 0.664115
- `C_delta` = 4.91

A formulação física da camada limite exige que a velocidade do escoamento tenda à velocidade do escoamento livre longe da placa. A simulação numérica atendeu a essa restrição com êxito: para um domínio suficientemente longo ($\eta_{max}$ = 20.0), o perfil adimensional de velocidades convergiu perfeitamente, registrando o valor alcançado de $f'(\eta_{max})$ = 1.000000. O erro do Método do Tiro para essa condição foi de 0.1161E-06, comprovando numericamente que $f'(\eta)$ vai para 1.

Com o valor da derivada segunda na parede, foi possível extrair a constante associada ao coeficiente de atrito local $Cf$ = 0.664115. A precisão dessa constante é fundamental do ponto de vista de engenharia, pois ela permite determinar o perfil real de tensão de cisalhamento na parede imposta por qualquer fluido Newtoniano, desde que as propriedades dimensionais do escoamento e o número de Reynolds local sejam conhecidos.

A espessura da camada limite resultou no coeficiente $C_{delta}$ = 4.91. Este resultado apresenta excelente concordância com o valor clássico da literatura ($C_{delta}$ = 4.92). A variação é reflexo da resolução da malha espacial, como a malha foi construída em intervalos fixos de $\delta\eta$ = 0.01, o nó exato que ocorre a transição de velocidade pode cair entre dois pontos discretizados, justificando a divergência posicional.

O método do tiro convergiu para o chute inicial da derivada segunda na parede no valor de $f"(0)$ = 0.332057. Este resultado replica com exatidão o valor clássico da literatura de referência, evidenciando a integridade e robustez do algoritmo. Apesar da precisão do resultado, vale destacar possíveis fontes de erro numérico inerentes à configuração escolhida para o solver:
* **Discretização e passo de integração:** o método RK4 tem erro de truncamento local de ordem $(\Delta\eta)^5$. Passos muito grandes subamostrariam a curva e impediriam a detecção precisa de variações da derivada. Passo muito pequenos geraram um excesso de iterações, abrindo margem para a degradação do resultado por conta do erro acumulado de arredondamento de precisão de ponto flutuante nas sucessivas somas de $\eta$.
* **Aproximação numérica do infinito ($\eta_{max}$):** adotar $\eta_{max}$ = 20.0 foi conservador e englobou o domínio necessário para o fluido estabilizar. Se o valor adotado fosse muito curto, o Método do Tiro seria forçado a induzir uma inclinação irreal na curvva para fazer o escoamento atingir a velocidade final prematuramente, o que corromperia o valor de $f"(0)$. Domínios excessivamente grandes apenas teriam um gasto computacional maior sem nenhum incremento de precisão física,pois as derivadas superiores se anulam ao longo de zonas distantes.

## Bibliografia
1. Roteiro PPC5 e APC5 disponibilizado pelo professor.
2. S. C. Chapra, R. P. Canale. “Métodos Numéricos para Engenharia.”, McGrawHill, 5a edição (2008): 1-825.