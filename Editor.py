# Editoração
from stringcolor import  *


def margem():
    print('-'*50)


def error(texto):
    print(cs(f'>>> ERRO : {texto} <<< ','white','red').bold(), end='\n')


def aviso(texto):
    print(cs(f'\n ATENÇÃO : {texto}','gold').bold(), end='\n')


def list(texto):
    print(cs(f'- {texto}', 'Cyan'), end='\n')


def endereco(texto):
    print(cs(f'\n-> {texto}','gold').bold(), end='\n')


def resposta(texto):
    print(cs(f' --- {texto}','grey2'), end='\n')


def titulo(texto,autor= 'Lucas Eduardo Silva Braga'):

    print('\n')
    print(bold('#'*50).cs('lime2'), end='\n')
    print(bold(f'{texto}').underline().cs('lime2'), end='\n')
    print(bold(f'Por : {autor}').cs('lime2'))
    print(bold('#' * 50).cs('lime2'), end='\n')

def relatorio_titulo(texto):
    print('\n')
    print(bold('-' * 50).cs('LightGrey10'), end='\n')
    print(bold(f'{texto.upper()}').cs('LightGrey10'), end='\n')
    print(bold('-' * 50).cs('LightGrey10'), end='\n')

def relatorio_item(texto):
    print(cs(f'- {texto}', 'LightGrey10'), end='\n')

def relatorio_end():
    print(bold('-' * 50).cs('LightGrey10'), end='\n')
    
def beep():
	import os

	duration = 0.05  #seg
	freq = 2000 #hz
	os.system(f'play -nq -t alsa synth {duration} sine {freq}')

