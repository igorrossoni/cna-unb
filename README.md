# Cálculo Numérico Aplicado (CNA) - Engenharia Mecânica

* **Instituição:** Universidade de Brasília (UnB)
* **Faculdade:** Faculdade de Tecnologia (FT)
* **Departamento:** Departamento de Engenharia Mecânica (ENM)
* **Estudante:** Igor de Oliveira Rossoni
* **Matrícula:** 222031279
* **Semestre Letivo:** 2026/1

Este é um repositório desenvolvido para a disciplina Cálculo Numérico Aplicado (CNA), buscando armazenar e organizar os algoritmos, relatórios e práticas desenvolvidas durante o semestre de 2026/1.

## Motivação para o estudo de métodos numéricos


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
├── PPC_03/                       # Diretório da Prática 02 (Método de Thomas)
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
├── PPC_04/                       # Diretório da Prática 02 (Maximização de 2 variáveis)
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
