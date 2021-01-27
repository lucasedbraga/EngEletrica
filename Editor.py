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
    print(bold(f'           {texto}          ').underline().cs('lime2'), end='\n')
    print(bold(f'            {autor}             ').underline().cs('lime2'))
    print(bold('#' * 50).cs('lime2'), end='\n')