#######################################
# Atividade 01 e 02 -  ENE103 2020/3  #
# Lucas Eduardo Silva Braga 201570023 #
#######################################

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

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
x = np.array(meses).flatten()
dados_velmed_ano = np.array(dados_velmed_ano)
legenda_ano = anos.astype(str)

##############################################################################################################

# 4- HISTOGRAMA COM A DISTRIBUIÇÃO ANUAL DA VELOCIDADE DO VENTO

#Iniciando o contador de horas/ano
contador = np.zeros((15))
#Coletando o número de Horas
for col in range(3, len(dados.columns)):
    for lin in range(len(dados)):
        velocidade = dados.iloc[lin][col]

        if velocidade < 0.5:
            contador[0]+= 1
        if 0.5 < velocidade < 1.5:
            contador[1] += 1
        if 1.5 < velocidade < 2.5:
            contador[2] += 1
        if 2.5 < velocidade < 3.5:
            contador[3] += 1
        if 3.5 < velocidade < 4.5:
            contador[4] += 1
        if 4.5 < velocidade < 5.5:
            contador[5] += 1
        if 5.5 < velocidade < 6.5:
            contador[6] += 1
        if 6.5 < velocidade < 7.5:
            contador[7] += 1
        if 7.5 < velocidade < 8.5:
            contador[8] += 1
        if 8.5 < velocidade < 9.5:
            contador[9] += 1
        if 9.5 < velocidade < 10.5:
            contador[10] += 1
        if 10.5 < velocidade < 11.5:
            contador[11] += 1
        if 11.5 < velocidade < 12.5:
            contador[12] += 1
        if 12.5 < velocidade < 13.5:
            contador[13] += 1
        if 13.5 < velocidade < 14.5:
            contador[14] += 1
        if 14.5 < velocidade < 15:
            contador[15] += 1

#Dividiando o número de horas por ano
contador = contador/(len(anos))

# PlOT
plt.figure()
plt.suptitle('HISTOGRAMAS - DADOS DE VENTO - MOSSORÓ (RN)')
plt.subplot(2, 1, 1)
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
plt.subplots_adjust(right=0.85, left=0.1)

plt.subplot(2, 1, 2)
plt.title('Distribuição Anual da Velocidade do Vento')
eixo_x = list(range(0,15))
plt.bar(eixo_x, contador, color='blue', label= 'Horas/Ano')
plt.legend(bbox_to_anchor=(1.0, 1.0), loc='best', borderaxespad=0.1)
plt.xlabel('Velocidade (m/s)')  # Legenda Eixo X
plt.ylabel('Horas/Ano')  # Legenda Eixo Y
plt.subplots_adjust(hspace=0.4, right=0.85, left=0.1)
plt.show()

##############################################################################################################




