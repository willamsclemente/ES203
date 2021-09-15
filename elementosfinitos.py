import numpy as np
import matplotlib.pyplot as plt

def matriz_local(): #Definindo a matriz local
    ke = np.zeros((2,2), dtype=np.int)
    ke[0][0] = -1
    ke[0][1] = 1
    ke[1][0] = 1
    ke[1][1] = -1
    #print ke # imprimindo a matriz local, quando necessario
    return ke

def matriz_global(N, ke): #Definindo a matriz global
    K = np.zeros((N,N), dtype=np.int)
    j = -1

    for i in range(0,N):
        if i == 0:
            K[i][0] = ke[0][0]
            K[i][1] = ke[0][1]
        elif i == N-1:
            K[i][N-2] = ke[1][0]
            K[i][N-1] = ke[1][1]
        else:
            K[i][j] = ke[1][0]
            K[i][j+1] = ke[1][1] + ke[0][0]
            K[i][j+2] = ke[0][1]
        
        j+=1
    #print K #Imprimindo a matriz global, quando necessario
    return K

def matriz_global_x_er(K, N1, N2, e0, er1, er2, L): #Multiplicando a matriz global pela permisividade absoluta
    for i in range(0, N1):
        for j in range(0, N1+N2):
            K[i][j] *= (e0*er1)/L
           
    
    for i in range(N1, N1+N2):
        for j in range(0, N1+N2):
            K[i][j] *= (e0*er2)/L
    
    key = 0
    for j in range(0, N1+N2):
        if K[N1][j] != 0 and key == 0:
            K[N1][j] = (e0*er1)/L
            K[N1][j+1] = -(e0*(er1+er2))/L
            K[N1][j+2] = (e0*er2)/L
            key = 1
    #print K #Imprimindo a matriz global multiplicada pela permisividade absoluta

def solucao_sistema(K, N1, N2, V0, VN_1): #Solucionando o sistema de equacoes e encontrando os potencias eletricos
    b = np.zeros((N1+N2), dtype=np.float64)
    b -= K[:,0]*V0
    K = np.delete(K, (0), axis=1)
    K = np.delete(K, (0), axis=0)
    b = np.delete(b, (0), axis=0)
   
    b -= K[:, N1+N2-2]*VN_1
    K = np.delete(K, (N1+N2-2), axis=1)
    K = np.delete(K, (N1+N2-2), axis=0)
    b = np.delete(b, (N1+N2-2), axis=0)
    
    V = np.linalg.solve(K,b)

    V = np.array([float(V0)]+list(V)+[float(VN_1)])
    #print V #Imprimindo os potencias eletricos encontrados, quando necessario 
    return V
    
def plot(V, N1, N2, d1, d2): #Plotagem do grafico dos potencias eletricos x distancia
    #inc = (float(d1 + d2) / float(N1 + N2 -1))
    inc1 = (float(d1)/float(N1))
    inc2 = (float(d2)/float(N2-1))
    z = []
   
    for i in range(0, N1+1):
        if i == N1:
            z += [d1]
        else:
            z  +=  [i*inc1]

    for i in range(1, N2):
        if i == N2-1:
                z += [d1+d2]
        else:
            z  +=  [d1 + (i*inc2)]

    z = np.array(z, dtype = float)
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.xlabel('Distancia(mm)')
    plt.ylabel('Potencial eletrico(V)')
    ax.plot(z,V, color = 'red')
    ax.set_title('V x z')
    plt.show()

def capacitancia(V, e0, er2, L, Lc, N): #Calculo e plotagem da capacitancia
    fig = plt.figure()
    ax = plt.subplot(111)
    A = Lc*Lc
    C = (e0*er2*A*(V[N-1]-V[N-2]))/(L*(V[N-1]-V[0]))
    #print C #Imprimindo capacitancia quando necessario
    plt.ylabel('Capacitancia(F)')
    ax.plot(0,C,'o',color = 'red')
    ax.set_title('Capacitancia')
    plt.show()


#Definindo os dados dos problemas e chamando as funcoes
e0 = 8.85*10**-12
er1 = 2
er2 = 4
N1 = 2
N2 = 3
N = N1+N2
d1 = 1.0
d2 = 1.0
L = ((d1 + d2)/(N-1))
Lc = 2*10**-2
V0 = 0
VN_1 = 1


ke = matriz_local()
K = matriz_global(N, ke)
K = K.astype(float)
matriz_global_x_er(K, N1, N2, e0, er1, er2, L)
V = solucao_sistema(K, N1, N2, V0, VN_1)
plot(V, N1, N2, d1, d2)
capacitancia(V, e0, er2, L, Lc, N)

