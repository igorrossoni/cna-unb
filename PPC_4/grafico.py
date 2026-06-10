# Geração do gráfico das curvas de nível e comparação dos dois métodos

import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# 1. Leitura dos Dados do Fortran
# =========================================================
# Lendo function.dat
dados_malha = np.loadtxt('function.dat')
X = dados_malha[:, 0]
Y = dados_malha[:, 1]
Z = dados_malha[:, 2]

# Dimensão original do grid (linhas e colunas)
pontos_y = len(np.unique(Y))
pontos_x = len(np.unique(X))

# Redimensionando os arrays 1D para matrizes 2D do grid
X_grid = X.reshape((pontos_x, pontos_y))
Y_grid = Y.reshape((pontos_x, pontos_y))
Z_grid = Z.reshape((pontos_x, pontos_y))

# Aclive Máximo (output1.dat)
# As colunas são: iter erro h x y dfx dfy (x é índice 3, y é índice 4)
am_data = np.loadtxt('output1.dat')
am_x = am_data[:, 3]
am_y = am_data[:, 4]

# Gradientes Conjugados (output2.dat)
gc_data = np.loadtxt('output2.dat')
gc_x = gc_data[:, 3]
gc_y = gc_data[:, 4]

# =========================================================
# 2. Configuração e Criação do Gráfico
# =========================================================
plt.figure(figsize=(10, 8))

# Desenhando as curvas de nível
# levels define a quantidade de "degraus" da montanha. Ajuste se quiser mais detalhes.
contorno = plt.contourf(X_grid, Y_grid, Z_grid, levels=25, cmap='viridis', alpha=0.8)
plt.colorbar(contorno, label='f(x,y)')

# Adicionando linhas pretas para destacar as curvas
plt.contour(X_grid, Y_grid, Z_grid, levels=25, colors='black', linewidths=0.5, alpha=0.5)

# Plotando a trajetória do Aclive Máximo
plt.plot(am_x, am_y, color='red', linestyle='--', marker='o', linewidth=2, markersize=5, label='Aclive Máximo')

# Plotando a trajetória dos Gradientes Conjugados
plt.plot(gc_x, gc_y, color='gray', linestyle='-.', marker='s', linewidth=2, markersize=6, label='Gradientes Conjugados')

# Marcando o ponto Ótimo Analítico com uma estrela
plt.plot(2, 1, color='gold', marker='*', markersize=15, markeredgecolor='black', label='Ótimo Analítico (2,1)')

# =========================================================
# 3. Estilização Final
# =========================================================
plt.title('Aclive Máximo vs Gradientes Conjugados', fontsize=14, fontweight='bold')
plt.xlabel('Eixo X', fontsize=12)
plt.ylabel('Eixo Y', fontsize=12)
plt.xlim([-2.5, 2.25])
plt.ylim([-1.0, 3.25])
plt.legend(loc='lower right', framealpha=0.9)
plt.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig('curvas_de_nivel.png', dpi=300, bbox_inches='tight')
