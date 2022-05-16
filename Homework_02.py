'''
Script done by:
Bruno Lopes 202000210
Gonçalo Cachado 202000190
Samuel Correia 202000094
Turma: Binf21
'''
from asyncore import write
from itertools import count
import sys
file_name = sys.argv[1]
file_outgroup = sys.argv[2]
file_ngen = sys.argv[3]

'''
Vai abrir e ler o ficheiro com o nome passado por argumento pela variável file_name.
As linhas do ficheiro vão ser guardados como itens da lista file_lines que
será retornada.
'''
def read_file(file_name):
    with open(file_name) as user_file:
        file_lines = user_file.readlines()
    user_file.close()
    return(file_lines)

'''
Vai receber a lista de linhas file_lines como argumento.
As linhas que contém o nome das sequências com mais de 99 caracteres
serão convertidos a nomes com 99 caracteres terminando sempre com um número.
Estas alterações serão guardadas na file_lines no lugar dos nomes anteriores e 
está variável será retornada.
'''
def name_shortener(file_lines):
    count1 = 0
    count2 = 0
    for line in file_lines:
        char_count = - line.count('\n')
        if((line[0] == '>') and ((len(line) + char_count)) > 99):
            line = line[0:100 - len(str(count2))] + str(count2) + '\n'
            file_lines[count1] = line
            count2 += 1
        count1 += 1
    return(file_lines)   

'''
Irá receber por argumento o file_lines e vai guardar os nomes das sequências
na variável seq_names, será também contado o numéro de sequências que existem e
guardado no num_seq.
O num_seq e seq_names serão retornados.
'''
def ntax_counter(file_lines):
    seq_names = [x for x in file_lines if x[0] == '>']
    num_seq = len(seq_names)
    return(num_seq,seq_names)

'''
Irá receber por argumento o file_lines e num_seq. 
Vai ser guardado as linhas das sequências e retirado os caracteres do file_lines \n
na variável seq_no_spaces. O num_char vai guardar o resultado da divisão do somatório de todos os caracteres das
sequências com o número de sequências o num_seq.
O num_char será retornado.
'''
def nchar_counter(file_lines,num_seq):    
    seq_no_spaces = [x.replace('\n','') for x in file_lines if x[0] != '>']
    num_char = int(sum([len(x) for x in seq_no_spaces])/num_seq)
    return(num_char)

'''
Vai receber o nome das sequências, seq_names por argumento.
Irá calcular o tamanho dos nomes do seq_names, que será guardado em
seq_names_length. É obtido o maior valor do seq_names_length e guardado no max_seq_length.
Os valores dos espaços a dar para alinhar as linhas será feito através da subtração do 
max_seq_length com os valores da seq_names_length, isto vai ser guardado na lista
seq_normalize_length que vai ser retornado. 

'''
def line_normalizer(seq_names):
    seq_names_length = [(len(x) - 2) for x in seq_names]
    max_seq_length = max(seq_names_length)
    seq_normalize_length = [(max_seq_length - x) for x in seq_names_length]
    return(seq_normalize_length)

'''
Irá receber as linhas do ficheiro, file_lines e a lista dos nomes, seq_names por argumento.
Será guardado no nex_format o resultado da transformação de fasta para nexus em string, para isso
vai ser iterado todas as linhas do file_lines e feitas alterações. Para alinhar as
novas linhas é usado a função line_normalizer e dado os nomes das sequências a alinhar.
As alterações de nexus para fasta serão dadas pelo retorno do nex_format.
'''
def nexus_transformer(file_lines,seq_names):
    nex_format = ''
    counter = 0
    for lines in file_lines:
        if(lines[0]== '>'):
            lines = ' ' * line_normalizer(seq_names)[counter] + lines
            lines = lines.replace('>',' ')
            lines = '\n' + lines.replace('\n','  ')
            nex_format = nex_format + lines 
            counter += 1
        else:
            nex_format = nex_format + lines.replace('\n','')
    return(nex_format)

'''
Será recebido por argumento as linhas de nexus através do nex_format, o ntax através do num_seq,
o nchar através do num_char, o outgroup por file_outgroup e ngen por file_ngen.
Será formato para nexus e guardado na variável nex_format que irá ser retornada.
'''
def nex_formatter(nex_format,num_seq,num_char,file_outgroup,file_ngen):
    nex_format = ('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX=' + str(num_seq) + ' NCHAR=' + str(num_char) + ';\n'
    + 'FORMAT DATATYPE=DNA MISSING=N GAP=-;\nMATRIX' + nex_format
    + '\n  ;\nEND;\n\nbegin mrbayes;\n  set autoclose=yes;'
    + '\n  outgroup ' + file_outgroup + '\n  mcmcp ngen=' +file_ngen + ' '
    + 'printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename=MyRun01;'
    + '\n  mcmc;\n  sumt filename=MyRun01;\nend;')
    return(nex_format)

file = read_file(file_name)
short_name = name_shortener(file)
ntax, seqs = ntax_counter(short_name)
nchar = nchar_counter(short_name,ntax)
nex = nexus_transformer(file,seqs)
print(nex_formatter(nex,ntax,nchar,file_outgroup,file_ngen))

