########################################################
# TVC 01 - Sistemas de Geração Éolica - ENE103 2020/3  #
# Lucas Eduardo Silva Braga - Matricula 201570023      #
########################################################
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
import Editor
from pathlib import Path
from scipy.interpolate import interp1d


############################################################################################################
Editor.titulo('CÁLCULO DA ENERGIA GERADA P/ MICROTURBINA EÓLICA  ')

Editor.aviso('Pressione ENTER Para começar')
input()

# 1- PLOTAR CURVA DE POTÊNCIA DA TURBINA: GERAR 246 V2
Editor.resposta('Plotando Curva de Potência da Microturbina: GERAR 246')
v_gerar_fornecido = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
velmax = max(v_gerar_fornecido)
pot_gerar_fornecido = [0, 8, 32, 72, 150, 240, 380, 500, 650, 805, 940, 1030, 1050, 980, 700]

plt.figure()
plt.title('Curva de Potência GERAR246')
plt.plot(v_gerar_fornecido, pot_gerar_fornecido)
plt.xlabel('Vel(m/s)')
plt.ylabel('Watt')
#plt.show()

# 2 - GRÁFICO INTERPOLADO
Editor.resposta('Plotando Curva Interpolada de Potência da Microturbina: GERAR 246')
x_interpol = np.linspace(0, velmax, num=10*velmax, endpoint=True)
pot_interpol = interp1d(v_gerar_fornecido, pot_gerar_fornecido, kind='cubic', fill_value='extrapolate')
plt.figure()
plt.title('Curva de Potência GERAR246')
plt.plot(v_gerar_fornecido, pot_gerar_fornecido)
plt.plot(x_interpol, pot_interpol(x_interpol), '--', color='red')
plt.legend(['Curva de Potência', 'Curva Interpolada'], loc='best')
plt.xlabel('Vel(m/s)')
plt.ylabel('Watt')
#plt.show()


# 3 - Curva de Potência Disponível

y_energia = []
x_energia = np.linspace(0, velmax, num=10*velmax)
ro = 1.225
passo = 0.1

for v in range(0, 10*velmax):
    v = passo*v
    Ed = (0.5*ro*(math.pi*((2.46/2)**2))*(v**3))
    y_energia.append(Ed)

betz = [(16/27)*y for y in y_energia]
pratica = [0.4*y for y in y_energia]



Editor.resposta(f'Plotando Curva de Energia Cinética Disponível ')
plt.figure()
plt.title('Curva de Energia Cinética Disponível')
plt.plot(x_energia, y_energia,  color='green', label='Cp = 100%')
plt.plot(x_energia, betz,  color='red', label='Cp = 59,26%')
plt.plot(x_energia, pratica,  color='yellow', label='Cp = 59,26%')
plt.plot(x_interpol, pot_interpol(x_interpol), '--', color='Blue')
plt.legend(['Potência Disponível - Cp=100%', 'Potência no Limite de Betz - Cp=59,26%','Curva Turbina Prática - Cp=40%','Potência da Turbina - 1000W'], loc='best')
plt.xlabel('Vel(m/s)')
plt.ylabel('Watt')
#plt.show()


# 4- LEITURA DO ARQUIVO DE DADOS


Editor.endereco('Escolha de Dados de Vento')
process = True
operacional = False
while process is True:
    pasta = Path('dados_vento/')
    print('\nLISTA DE CIDADES: ')
    for arquivos in pasta.iterdir():
        Editor.list(f'{arquivos.name}')

    Editor.aviso('Para sair digite CANCEL')
    sistema = input('\n Digite o nome do sistema: ->')

    if sistema.upper() == 'CANCEL':
        process = False
        break

    try:
        # Abrindo Arquivo de Dados
        try:
            arquivo_csv = pd.read_csv(f'dados_vento/{sistema}.csv', sep='\t', header=1)
            # Criando Banco de Dados
            dados = []
            for linha in range(len(arquivo_csv)):

                amostra = arquivo_csv.iloc[linha, :26]
                leitura_data = amostra['HORA UTC']

                leitura_dia = leitura_data[0:2]
                leitura_mes = leitura_data[3:6]
                leitura_ano = leitura_data[7:11]

                leitura = amostra[1:]
                vel_alt_exp = leitura
                vel_alt_log = leitura
                for pos in range(len(leitura)):
                    elemento = leitura[pos]
                    if type(elemento) == str:
                        valor = elemento.replace(',', '.')
                        elemento = float(valor)
                    leitura[pos] = elemento


                leitura = np.array(leitura)
                vel_alt_exp = np.array(vel_alt_exp)
                vel_alt_log = np.array(vel_alt_log)
                data = np.array([leitura_dia, leitura_mes, leitura_ano])
                dado = np.concatenate((data, leitura))

                dados.append(dado)

            dados = pd.DataFrame(dados)

            Editor.resposta(f'Carregando arquivo {sistema}.csv')
############################################################################################################

            # 5- TRATAMENTO DE DADOS FALTANTES

            Editor.resposta('Tratando dados faltantes')
            check_of_nan = dados.isnull()
            print('\n')
            for col in range(3, len(dados.columns)):
                for lin in range(len(dados)):
                    if math.isnan(dados[col][lin]):
                        # Percorrer o valores de Mês/Dia
                        flag_mes = dados.iloc[lin, 1]
                        flag_dia = dados.iloc[lin, 0]
                        flags = []
                        dados_calculo = []
                        for check in range(0, len(dados)):
                            check_mes = dados.iloc[check, 1]
                            check_dia = dados.iloc[check, 0]
                            if check_mes == flag_mes:
                                if check_dia == flag_dia:
                                    flags.append(check)

                        for flag in flags:
                            dados_calculo.append(dados.iloc[flag, col])
                        dados_calculo = pd.DataFrame(dados_calculo)
                        # Calculo da Média de Velocidade na Data
                        media = dados_calculo[0].mean()
                        # Substituição
                        dados[col][lin] = media

            process = False
            operacional = True
        except:
            Editor.error('Arquivo Inválido')
    except:
        Editor.aviso('Comando não reconhecido, tente novamente!')




###############################################################################
# CALCULO DAS ALTURAS
exp = ((50/10)**0.20)
print(exp)
log = ((math.log(50/0.05))/(math.log(10/0.05)))
dados_alt_exp = exp * np.array(dados.iloc[:, 3:], dtype=float)
dados_alt_log = log * np.array(dados.iloc[:, 3:], dtype=float)
dados_alt_exp = pd.DataFrame(dados_alt_exp)
dados_alt_log = pd.DataFrame(dados_alt_log)
data = dados.iloc[:,:3]
dados_alt_exp = pd.concat((data, dados_alt_exp), axis=1, ignore_index=True)
dados_alt_log = pd.concat((data, dados_alt_log), axis=1, ignore_index=True)
#dados_alt_exp = pd.DataFrame(dados_alt_exp)
#dados_alt_log = pd.DataFrame(dados_alt_log)

print(dados_alt_log)
print(dados_alt_exp)

###############################################################################




# 4- Agrupando os dados
if operacional:
    anos = dados.iloc[:][2]
    anos = pd.DataFrame(anos).drop_duplicates()
    meses = dados.iloc[:][1]
    meses = pd.DataFrame(meses).drop_duplicates()
    anos = anos.reset_index(drop=True)
    meses = meses.reset_index(drop=True)
    anos = anos.astype(str)
    meses = meses.astype(str)
    legenda_ano = anos.astype(str)

    # 5 - Cálculos Estatísticos

    dados_velmed_ano = []
    dados_desvpad_ano = []
    dados_velmed_ano_exp = []
    dados_desvpad_ano_exp = []
    dados_velmed_ano_log = []
    dados_desvpad_ano_log = []
    for ano in range(len(anos)):
        dados_velmed = []
        dados_desvpad = []
        dados_velmed_exp = []
        dados_desvpad_exp = []
        dados_velmed_log = []
        dados_desvpad_log = []
        for mes in range(len(meses)):
            #Selecionando conjunto de dados com Mês/Ano
            dado_med = dados.loc[dados[2] == anos.iloc[ano][2]]
            dado_med = dado_med.loc[dado_med[1] == meses.iloc[mes][1]]
            dado_med_exp = dados_alt_exp.loc[dados_alt_exp[2] == anos.iloc[ano][2]]
            dado_med_exp = dado_med_exp.loc[dado_med_exp[1] == meses.iloc[mes][1]]
            dado_med_log = dados_alt_log.loc[dados_alt_log[2] == anos.iloc[ano][2]]
            dado_med_log = dado_med_log.loc[dado_med_log[1] == meses.iloc[mes][1]]
            med = dado_med.iloc[:, 3:]
            med = pd.DataFrame(med)
            med_exp = dado_med_exp.iloc[:, 3:]
            med_exp = pd.DataFrame(med_exp)
            med_log = dado_med_log.iloc[:, 3:]
            med_log = pd.DataFrame(med_log)
            #Cálculo da Média de Velocidade por Hora
            med = med.mean()
            med_exp = med_exp.mean()
            med_log = med_log.mean()
            #Cálculo de Desvio Padrão do Mês
            desvpad = med.std()
            #Cálculo da Média de Velocidade por Mês
            vel_med = med.mean()
            vel_med_exp = med_exp.mean()
            vel_med_log = med_log.mean()
            #Agrupando os Dados do Mês
            dados_velmed.append(vel_med)
            dados_velmed_exp.append(vel_med_exp)
            dados_velmed_log.append(vel_med_log)
            dados_desvpad.append(desvpad)
        #Agrupando os Dados do Ano
        dados_velmed_ano.append(dados_velmed)
        dados_velmed_ano_exp.append(dados_velmed_exp)
        dados_velmed_ano_log.append(dados_velmed_log)
        #Desvio Padrão do Ano
        dados_desvpad_ano.append(dados_desvpad)

    dados_velmed_ano = pd.DataFrame(dados_velmed_ano)
    dados_velmed_ano_exp = pd.DataFrame(dados_velmed_ano_exp)
    dados_velmed_ano_log = pd.DataFrame(dados_velmed_ano_log)
    dados_desvpad_ano = pd.DataFrame(dados_desvpad_ano)
    velmed_ano = []
    velmed_ano_exp = []
    velmed_ano_log = []
    desvpad_ano = []

    for ano in range(len(dados_velmed_ano)):
        velmed = dados_velmed_ano.iloc[ano, :].mean()
        velmed_exp = dados_velmed_ano_exp.iloc[ano, :].mean()
        velmed_log = dados_velmed_ano_log.iloc[ano, :].mean()
        desvpad = dados_desvpad_ano.iloc[ano, :].mean()
        desvpad_ano.append(desvpad)
        velmed_ano.append(velmed)
        velmed_ano_exp.append(velmed_exp)
        velmed_ano_log.append(velmed_log)
    x = np.array(meses).flatten()
    dados_velmed_ano = np.array(dados_velmed_ano)
    dados_velmed_ano_exp = np.array(dados_velmed_ano_exp)
    dados_velmed_ano_log = np.array(dados_velmed_ano_log)




    # Distribuição de Rayleigh
    c_list = []
    c_list_exp = []
    c_list_log = []

    for ano in range(len(velmed_ano)):
        c = (2/math.sqrt(math.pi))*velmed_ano[ano]
        c_list.append(c)
        c_exp = (2 / math.sqrt(math.pi)) * velmed_ano_exp[ano]
        c_list_exp.append(c_exp)

        c_log = (2 / math.sqrt(math.pi)) * velmed_ano_log[ano]
        c_list_log.append(c_log)




    f_rayleigh = []
    energia_disp_ano = []
    energia_gerada_ano = []
    energia_cap_ano = []

    f_rayleigh_exp = []
    energia_disp_ano_exp = []
    energia_gerada_ano_exp = []
    energia_cap_ano_exp = []

    f_rayleigh_log = []
    energia_disp_ano_log = []
    energia_gerada_ano_log = []
    energia_cap_ano_log = []

    eixo_x_rey = np.arange(0.0, velmax, passo)

    for ano in range(len(c_list)):
        f_ano = []
        energia_disp = 0
        energia_cap = 0
        energia_ger = 0
        count = 0

        f_ano_exp = []
        energia_disp_exp = 0
        energia_cap_exp = 0
        energia_ger_exp = 0
        count_exp = 0

        f_ano_log = []
        energia_disp_log = 0
        energia_cap_log = 0
        energia_ger_log = 0
        count_log = 0


        for v in range(0, 10*velmax):
            v = passo*v
            rayleigh = ((2*v) / ((c_list[ano])**2)) * (math.e**(-1*((v/c_list[ano])**2)))/10
            energia_disp += (0.5 * ro * (math.pi * ((2.46 * 0.5) ** 2)) * (v ** 3))*(rayleigh*8760)
            energia_ger += pot_interpol(v)*(rayleigh*8760)
            count += rayleigh

            f_ano.append(rayleigh)

            rayleigh_exp = ((2 * v) / ((c_list_exp[ano]) ** 2)) * (math.e ** (-1 * ((v / c_list_exp[ano]) ** 2))) / 10
            energia_disp_exp += (0.5 * ro * (math.pi * ((2.46 * 0.5) ** 2)) * (v ** 3)) * (rayleigh_exp * 8760)
            energia_ger_exp += pot_interpol(v) * (rayleigh_exp * 8760)
            count_exp += rayleigh_exp

            f_ano_exp.append(rayleigh_exp)

            rayleigh_log = ((2 * v) / ((c_list_log[ano]) ** 2)) * (math.e ** (-1 * ((v / c_list_log[ano]) ** 2))) / 10
            energia_disp_log += (0.5 * ro * (math.pi * ((2.46 * 0.5) ** 2)) * (v ** 3)) * (rayleigh_log * 8760)
            energia_ger_log += pot_interpol(v) * (rayleigh_log * 8760)
            count_log += rayleigh_log

            f_ano.append(rayleigh_log)

        print(count)
        print(count_exp)
        print(count_log)

        energia_cap_ano.append(energia_cap)
        energia_gerada_ano.append(energia_ger)
        energia_disp_ano.append(energia_disp)
        f_rayleigh.append(f_ano)

        energia_cap_ano_exp.append(energia_cap_exp)
        energia_gerada_ano_exp.append(energia_ger_exp)
        energia_disp_ano_exp.append(energia_disp_exp)
        f_rayleigh_exp.append(f_ano_exp)

        energia_cap_ano_log.append(energia_cap_log)
        energia_gerada_ano_log.append(energia_ger_log)
        energia_disp_ano_log.append(energia_disp_log)
        f_rayleigh_log.append(f_ano_log)


    for ano in range(len(anos)):

        energia_disponivel = ((energia_disp_ano[ano]) / 10**6)
        energia_gerada = ((energia_gerada_ano[ano])/10**6)
        energia_capaz = ((pot_interpol(12.5)*8760)/10**6)
        fator_capacidade = (energia_gerada/energia_capaz)*100
        rendimento = (energia_gerada/energia_disponivel)*100

        ganho_energia_disponivel_exp = (energia_disp_ano_exp[ano]/energia_disp_ano[ano])*100
        ganho_energia_gerada_exp = (energia_gerada_ano_exp[ano]/energia_gerada_ano[ano])*100

        ganho_energia_disponivel_log = (energia_disp_ano_log[ano] / energia_disp_ano[ano])*100
        ganho_energia_gerada_log = (energia_gerada_ano_log[ano] / energia_gerada_ano[ano])*100


        ###########################################################################################################
        Editor.relatorio_titulo(f' RELATÓRIO - ENERSUD GERAR 246 ({sistema.upper()}: {legenda_ano.iloc[ano][2]})')
        Editor.relatorio_item(f'Energia Cinética Disponível = {energia_disponivel:.3f} MWh')
        Editor.relatorio_item(f'Energia Produzida Durante o Ano = {energia_gerada:.3f} MWh')
        Editor.relatorio_item(f'Fator de Capacidade = {fator_capacidade:.2f} %')
        Editor.relatorio_item(f'Relação Geração/Disponibilidade = {rendimento:.2f} %')
        Editor.relatorio_end()
        ###########################################################################################################
        Editor.relatorio_titulo(f'  CORREÇÃO DA ALTURA 50M PARA O CASO - RURAL')
        Editor.relatorio_item(f'Ganho Energia Cinética Disp- Exponencial = {ganho_energia_disponivel_exp:.2f} %')
        Editor.relatorio_item(f'Ganho Energia Produzida- Exponencial = {ganho_energia_gerada_exp:.2f} %')
        Editor.relatorio_item(f'Ganho Energia Cinética Disp- Logaritmo = {ganho_energia_disponivel_log:.2f} %')
        Editor.relatorio_item(f'Ganho Energia Produzida- Logaritmo = {ganho_energia_gerada_log:.2f} %')
        Editor.relatorio_end()


Editor.resposta('FIM')
##############################################################################################################
