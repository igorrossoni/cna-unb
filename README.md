# Cálculo Numérico Aplicado (CNA) - Engenharia Mecânica

* **Instituição:** Universidade de Brasília (UnB)
* **Faculdade:** Faculdade de Tecnologia (FT)
* **Departamento:** Departamento de Engenharia Mecânica (ENM)
* **Estudante:** Igor de Oliveira Rossoni
* **Matrícula:** 222031279
* **Semestre Letivo:** 2026/1

Este é um repositório desenvolvido para a disciplina Cálculo Numérico Aplicado (CNA), buscando armazenar e organizar os algoritmos, relatórios e práticas desenvolvidas durante o semestre de 2026/1.

## Motivação para o estudo de métodos numéricos
A engenharia moderna lida constantemente com desafios de alta complexidade que, na grande maioria das vezes, não podem ser solucionados apenas por vias analíticas tradicionais. O estudo de Métodos Numéricos surge como a ponte fundamental entre os modelos matemáticos teóricos e a aplicação prática, permitindo a resolução de sistemas reais — desde o dimensionamento de eixos de transmissão e análise de integridade estrutural até simulações avançadas em termodinâmica e mecânica dos fluidos.

A complexidade dos problemas contemporâneos exige o uso de métodos computacionais robustos e a colaboração estreita entre profissionais. Nesse cenário, a programação deixa de ser apenas uma habilidade complementar e torna-se uma ferramenta essencial e um requisito operacional para a resolução de problemas em Cálculo Numérico. O domínio dessas técnicas e o entendimento da lógica numérica "sob o capô" permitem modelar fenômenos com alta precisão, evitando que erros ocultos em softwares proprietários (caixas-pretas) comprometam cálculos estruturais ou simulações críticas.

Além disso, a adoção de boas práticas de estruturação e a filosofia de código aberto (Open Source) são vitais para a reprodutibilidade no cálculo científico. Ao desenvolver e compartilhar rotinas numéricas claras e validadas, o conhecimento deixa de ser um ativo privado e torna-se um bem coletivo. Isso garante que os algoritmos possam ser verificados e aprimorados globalmente, promovendo a excelência técnica, o rigor matemático e a inovação na engenharia.

## Estrutura do repositório
```text
cna-unb/
├── README.md                     # Metadados, motivação e bibliografia geral
├── livro-cna.pdf                 # Livro utilizado de referência para o estudo dos métodos numéricos
├── PPC_01/                       # Diretório da Prática 01 (Runge-Kutta de 4° ordem)
│   ├── README.md                 # Resumo, dicionário de variáveis, I/O e validação
│   ├── PPC1 e APC1.pdf           # Roteiro da atividade
│   ├── APC_1_igor_rossoni.pdf    # APC finalizada
│   ├── ppc1-rk4.f90              # Script principal com o algoritmo numérico
│   ├── plot1.gnu                 # Script de plotagem dos gráficos
│   └── outputs/                  # Saídas numéricas, logs e gráficos
|
├── PPC_02/                       # Diretório da Prática 02 (Método de Bairstow)
│   ├── README.md                 # Resumo, dicionário de variáveis, I/O e validação
│   ├── PPC2 e APC2.pdf           # Roteiro da atividade
│   ├── APC_2_igor_rossoni.pdf    # APC finalizada
│   ├── ppc2-bairstow.f90         # Script principal com o algoritmo numérico
│   ├── varredura.f90             # Script de varredura no plano bidimensional
│   ├── gerar_fractal.f90         # Script de plotagem do fractal
│   └── outputs/                  # Saídas numéricas, logs e gráficos
│       └── dados_varredura.csv
│       └── fractal_bairstow.png
|
├── PPC_03/                       # Diretório da Prática 03 (Método de Thomas)
│   ├── README.md                 # Resumo, dicionário de variáveis, I/O e validação
│   ├── PPC3 e APC3.pdf           # Roteiro da atividade
│   ├── APC_3_igor_rossoni.pdf    # APC finalizada
│   ├── ppc3-thomas.f90           # Script principal com o algoritmo numérico
│   ├── campo_com_geracao.gnu     # Script de plotagem do campo com geração interna
│   ├── evolucao_com_geracao.gnu  # Script de plotagem da evolução da temperatura com o tempo com geração interna
│   ├── campo_sem_geracao.gnu     # Script de plotagem do campo sem geração interna
│   ├── evolucao_sem_geracao.gnu  # Script de plotagem da evolução da temperatura com o tempo sem geração interna
│   └── outputs/                  # Saídas numéricas, logs e gráficos
│       └── campo_com_geracao.png
│       └── evolucao_com_geracao.png
│       └── campo_sem_geracao.png
│       └── evolucao_sem_geracao.png
│       └── saida.dat
|
├── PPC_04/                       # Diretório da Prática 04 (Maximização de 2 variáveis)
│   ├── README.md                 # Resumo, dicionário de variáveis, I/O e validação
│   ├── PPC4 e APC4.pdf           # Roteiro da atividade
│   ├── APC_4_igor_rossoni.pdf    # APC finalizada
│   ├── ppc4.f90                  # Script principal com o algoritmo numérico
│   ├── grafico.py                # Script de plotagem do gráfico em python
│   └── outputs/                  # Saídas numéricas, logs e gráficos
│       └── curvas_de_nivel.png
│       └── function.dat
│       └── output1.dat
│       └── output2.dat
|
├── PPC_05/                       # Diretório da Prática 05 (Equação de Blasius)
│   ├── README.md                 # Resumo, dicionário de variáveis, I/O e validação
│   ├── PPC5 e APC5.pdf           # Roteiro da atividade
│   ├── APC_5_igor_rossoni.pdf    # APC finalizada
│   ├── blasius.f90               # Script principal com o algoritmo numérico
│   ├── plot.gnu                  # Script de plotagem do gráfico Gnuplot
│   └── outputs/                  # Saídas numéricas, logs e gráficos
│       └── perfis.png
│       └── saida.dat
|
├── PPC_06/                       # Diretório da Prática 06 (Transferência de calor em aleta)
│   ├── README.md                 # Resumo, dicionário de variáveis, I/O e validação
│   ├── PPC6 e APC6.pdf           # Roteiro da atividade
│   ├── APC_6_igor_rossoni.pdf    # APC finalizada
│   ├── ppc6.f90                  # Script principal com o algoritmo numérico
│   ├── plot.gnu                  # Script de plotagem do gráfico Gnuplot
│   └── outputs/                  # Saídas numéricas, logs e gráficos
│       └── erro_percentual.png
│       └── linha_central.png
│       └── mapa_temp.png
│       └── centro.dat
│       └── gauss.dat
│       └── liebmann.dat
│       └── relax.dat
|
└── trabalhos_sala/           # Projetos feitos em sala
    ├── README.md             # Resumo do que tem no diretório e dos métodos
    ├── NewtonRaphson/
    ├── bisseccao-falsaPosicao/
    ├── decomp-LU-inversa/
    ├── eliminacao/           # Eliminação Gaussina
    ├── euler-heun/
    ├── interp-quad-razao-aurea/
    ├── liebmann/
    ├── muller
    ├── tiro/
    └── plot-graficos-python/ # Programas de plotagem de alguns gráficos em python
```

## Referências
* Métodos numéricos para Engenharia, Steven C. Chapra e Raymond P. Canale, McGrawHill;
