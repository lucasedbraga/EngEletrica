#######################################
# Atividade 03 e 04 -  ENE103 2020/3  #
# Lucas Eduardo Silva Braga 201570023 #
#######################################

import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
import statistics

############################################################################################################

# 1- LEITURA DO ARQUIVO DE DADOS

# Abrindo Arquivo de Dados
arquivo_csv = pd.read_csv('ventos_mossoro.csv', sep='\t', header=1)
# Criando Banco de Dados
dados = []
for linha in range(len(arquivo_csv)):

    amostra = arquivo_csv.iloc[linha, :]

    leitura_data = amostra['HORA UTC']
    leitura_dia = leitura_data[0:2]
    leitura_mes = leitura_data[3:6]
    leitura_ano = leitura_data[7:11]

    leitura = amostra[1:]
    for pos in range(len(leitura)):
        elemento = leitura[pos]
        if type(elemento) == str:
            valor = elemento.replace(',', '.')
            elemento = float(valor)
        leitura[pos] = elemento
    leitura = np.array(leitura)
    data = np.array([leitura_dia,leitura_mes,leitura_ano])
    dado = np.concatenate((data,leitura))
    dados.append(dado)

dados = pd.DataFrame(dados)
print(f' DADOS OBTIDOS PELA TABELA EXCEL: \n \n {dados}')
##############################################################################################################

# 2- TRATAMENTO DE DADOS FALTANTES

#Encontrando valores NaN
check_of_nan = dados.isnull()
for col in range(3, len(dados.columns)):
    for lin in range(len(dados)):
        if math.isnan(dados[col][lin]):
            # Percorrer o valores de Mês/Dia
            flag_mes = dados.iloc[lin,1]
            flag_dia = dados.iloc[lin,0]
            flags = []
            dados_calculo = []
            for check in range(0, len(dados)):
                check_mes = dados.iloc[check,1]
                check_dia = dados.iloc[check,0]
                if check_mes == flag_mes:
                    if check_dia == flag_dia:
                        flags.append(check)

            for flag in flags:
                dados_calculo.append(dados.iloc[flag,col])
            dados_calculo= pd.DataFrame(dados_calculo)
            #Calculo da Média de Velocidade na Data
            media = dados_calculo[0].mean(skipna=True)
            #Substituição
            dados[col][lin] = media
##############################################################################################################

# 3- PLOTAR HISTOGRAMAS COM A MÉDIA MENSAL

print(f' DADOS APÓS O TRATAMENTO: \n \n {dados}')

anos = dados.iloc[:][2]
anos = pd.DataFrame(anos).drop_duplicates()
meses = dados.iloc[:][1]
meses = pd.DataFrame(meses).drop_duplicates()

anos = anos.reset_index(drop=True)
meses = meses.reset_index(drop=True)

anos = anos.astype(str)
meses = meses.astype(str)

dados_velmed_ano = []
for ano in range(len(anos)):
    dados_velmed = []
    for mes in range(len(meses)):
        #Selecionando conjunto de dados com Mês/Ano
        dado_med = dados.loc[dados[2] == anos.iloc[ano][2]]
        dado_med = dado_med.loc[dado_med[1] == meses.iloc[mes][1]]
        med = dado_med.iloc[:, 3:]
        med = pd.DataFrame(med)
        #Cálculo da Média de Velocidade por Hora
        med = med.mean()
        #Cálculo da Média de Velocidade por Mês
        vel_med = med.mean()
        #Agrupando os Dados do Mês
        dados_velmed.append(vel_med)
    #Agrupando os Dados do Ano
    dados_velmed_ano.append(dados_velmed)


#Peparando dados para serem plotados
dados_velmed_ano = pd.DataFrame(dados_velmed_ano)
velmed_ano = []
desvpad_ano = []
for ano in range(len(dados_velmed_ano)):
    velmed = dados_velmed_ano.iloc[ano,:].mean()
    desvpad = statistics.stdev(dados_velmed_ano.iloc[ano,:])
    desvpad_ano.append(desvpad)
    velmed_ano.append(velmed)
x = np.array(meses).flatten()
dados_velmed_ano = np.array(dados_velmed_ano)
legenda_ano = anos.astype(str)

##############################################################################################################

# 4- HISTOGRAMA COM A DISTRIBUIÇÃO ANUAL DA VELOCIDADE DO VENTO



dados_hist = []
for ano in range(len(anos)):
    dados_hist_ano = []
    #Selecionando conjunto de dados com Ano
    dado_hist = dados.loc[dados[2] == anos.iloc[ano][2]]
    dado_hist = dado_hist.iloc[:, 3:]
    dado_hist = pd.DataFrame(dado_hist)
    dados_hist.append(dado_hist)


#Iniciando o contador de horas/ano
contador_hist = []
passo = 0.5
eixo_x = np.arange(0.0, 15.0, passo)

for ano in range(len(dados_hist)):
    contador = np.zeros((30))
    #Coletando o número de Horas
    for col in range(3, len(dados_hist[ano].columns)):
        for lin in range(len(dados_hist[ano])):
            for pos in range(len(contador)):
                velocidade = dados_hist[ano].iloc[lin][col]
                if int((0.5*(pos-passo)).__abs__()) < velocidade <= 0.5*(pos +0.5) :
                    contador[pos] += 1
                    #Dividiando o número de horas por ano
    contador_hist.append(contador)
contador_hist = pd.DataFrame(contador_hist)



# Distribuição de Rayleigh
c_list = []
for ano in range(len(velmed_ano)):
    c = (2/math.sqrt(math.pi))*velmed_ano[ano]
    c_list.append(c)


f_rayleigh = []
passo_dist = 0.1
eixo_x_rey = np.arange(0.0, 15.0, passo_dist)
for ano in range(len(c_list)):
    f_ano = []
    for v in range(0,150):
        v = passo_dist*v
        rayleigh = ((2*v)/((c_list[ano])**2))*(math.e**(-(v/c_list[ano])**2))
        f_ano.append(rayleigh)
    f_rayleigh.append(f_ano)



# Distribuição de Weibull

f_weibull = []

eixo_x_weib = np.arange(0.0, 15.0, passo_dist)
for ano in range(len(c_list)):
    f_ano = []
    for v in range(0,150):
        v = passo_dist*v
        k = (desvpad_ano[ano]/velmed_ano[ano])**(-1.086)
        c = (velmed_ano[ano])/(math.gamma(1+(1/k)))
        weibull = (k/c) * ((v/c)**(k-1)) * (math.e**(-((v/c)**k)))
        f_ano.append(weibull)

    f_weibull.append(f_ano)




# PlOT - Todos os Meses
plt.figure()
plt.suptitle('DADOS DE VENTO - MOSSORÓ (RN) : Comparação entre os Anos')
plt.title('Velocidade Média Mensal')
plt.plot(x, dados_velmed_ano[0], color='red', label=legenda_ano.iloc[0][2])
plt.plot(x, dados_velmed_ano[1], color='blue', label=legenda_ano.iloc[1][2])
plt.plot(x, dados_velmed_ano[2], color='yellow', label=legenda_ano.iloc[2][2])
plt.plot(x, dados_velmed_ano[3], color='green', label=legenda_ano.iloc[3][2])
plt.plot(x, dados_velmed_ano[4], color='#00BFFF', label=legenda_ano.iloc[4][2])
plt.plot(x, dados_velmed_ano[5], color='#00FF00', label=legenda_ano.iloc[5][2])
plt.legend(bbox_to_anchor=(1.0, 1.0), loc='best', borderaxespad=0.1)
plt.xlabel('Meses')  # Legenda Eixo X
plt.ylabel('Velocidade Média (m/s)')  # Legenda Eixo Y
plt.show()


for ano in range(len(anos)):
    # PlOT - Detalhado
    plt.suptitle(f'HISTOGRAMAS DE VENTO DETALHADOS - MOSSORÓ (RN) - ANO : {legenda_ano.iloc[ano][2]}')
    plt.subplot(2, 2, 1)
    plt.title('Velocidade Média Mensal')
    plt.plot(x, dados_velmed_ano[ano], color='red', label=legenda_ano.iloc[ano][2])
    plt.legend(bbox_to_anchor=(1.0, 1.0), loc='best', borderaxespad=0.1)
    plt.xlabel('Meses')  # Legenda Eixo X
    plt.ylabel('Velocidade Média (m/s)')  # Legenda Eixo Y
    plt.subplots_adjust(hspace=0.5, right=0.95, left=0.04)

    ax = plt.subplot(2, 2, 2)
    plt.title('Distribuição Anual da Velocidade do Vento')
    histogramas = ax.bar(eixo_x, contador_hist.iloc[ano,:], width=(passo), color='blue', label= legenda_ano.iloc[ano][2])
    for bar in histogramas:
        bar.set_edgecolor("black")
        bar.set_linewidth(0.5)
    plt.legend(bbox_to_anchor=(1.0, 1.0), loc='best', borderaxespad=0.1)
    plt.xlabel('Velocidade (m/s)')  # Legenda Eixo X
    plt.ylabel('Horas/Ano')  # Legenda Eixo Y
    plt.subplots_adjust(hspace=0.5, right=0.98, left=0.04)


    plt.subplot(2, 2, 3)
    plt.title(f'Distribuição de Rayleigh')
    plt.bar(eixo_x_rey, f_rayleigh[ano], width=(passo_dist/2), color='green', label=legenda_ano.iloc[ano][2])
    plt.legend(bbox_to_anchor=(1.0, 1.0), loc='best', borderaxespad=0.1)
    plt.xlabel('Velocidade (m/s)')  # Legenda Eixo X
    plt.ylabel('Probabilidade')  # Legenda Eixo Y
    plt.subplots_adjust(hspace=0.5, right=0.95, left=0.04)


    plt.subplot(2, 2, 4)
    plt.title(f'Distribuição de Weibull')
    plt.bar(eixo_x_weib, f_weibull[ano], width=(passo_dist / 2), color='#00BFFF', label=legenda_ano.iloc[ano][2])
    plt.legend(bbox_to_anchor=(1.0, 1.0), loc='best', borderaxespad=0.1)
    plt.xlabel('Velocidade (m/s)')  # Legenda Eixo X
    plt.ylabel('Probabilidade')  # Legenda Eixo Y
    plt.subplots_adjust(hspace=0.5, right=0.95, left=0.04)

    plt.show()

##############################################################################################################




