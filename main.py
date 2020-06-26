#!/usr/bin/env python
# -*- coding: UTF-8 -*-

##########Módulos##########

import sys
from Bio import Seq
from Bio import SeqIO
import re
import os
import blast
import muscle
import prosite_mod
from datetime import datetime
from datetime import date
import shutil

#########Funciones#########

def ayudalargo():

        '''
        Muestra un mensaje de ayuda detallado.
        '''

        print("\n\033[1;36m"+"DESCRIPCIÓN:"+"\033[00m"+"\033[;36m"
                +" Este script realiza las siguientes"+" tareas:\n")
        print("\t- BlastP de cada query aportada por el usuario frente"
                +" a proteínas obtenidas de Genbank.")
        print("\t- Alineamiento de los hits procedentes del BlastP"
                +" (usando Muscle).")
        print("\t- Árbol filogenético por el metodo Neighbor-Joining para"
                +" cada proteína query (usando Muscle).")
        print("\t- Búsqueda de dominios protéicos presentes en database de"
                +" Prosite indicando sus\n\t  características principales (name," 
                +" accession, description, domain...).")
        print("\n\033[1;36m"+"USO: "+"\033[00m"+"\033[;36m"
                +"./main.py [query] [dir_subject/] [identity cut-off]"
                +" [coverage cut-off] [prosite database]")
        print("\033[1;36m_________________________________________\033[00m")
        print("\tArgumentos a escribir:")
        print("\tEl primero es el query a analizar. Debe tener formato fasta.")
        print("\tEl segundo es el dir_subject, el directorio donde se encuentran"
                +" los subject,\n\tlos cuales deben tener formato Genbank.")
        print("\tEl tercero es el identity cut-off. Debe ser un número"
                +" entre 0 y 100.")
        print("\tEl cuarto es el coverage cut-off. Debe ser un número"
                +" entre 0 y 100.")
        print("\tEl quinto es la base de datos de PROSITE. Debe tener"
                +" formato Prosite.")        
        print("\033[1;36m_________________________________________\033[00m")


def ayudacorto ():

        """
        Muestra un mensaje de ayuda corto.
        """

        print("\n\033[1;36m"+"USO: "+"\033[00m"+"\033[;36m"
                +"./main.py [query] [dir_subject/] [identity cut-off]"
                +" [coverage cut-off] [prosite database]")
        print("\033[1;36m_________________________________________\033[00m")
        print("\tArgumentos:")
        print("\t[query]: secuencias a analizar. Debe tener formato fasta.")
        print("\t[dir_subject/]: directorio donde se encuentran los subject,"
                +" los cuales deben tener formato Genbank.")
        print("\t[identity cut-off]: número entre 0 y 100.")
        print("\t[coverage cut-off]: número entre 0 y 100.")
        print("\t[prosite database]: base de datos de PROSITE."
                +" Debe tener formato Prosite.")  
        print("\033[1;36m_________________________________________\033[00m")


def ayuda():

        '''
        '''

        for t in sys.argv:
                if t == "-h" or t == "--help":
                        ayudalargo()   
                        sys.exit()



def error():

        '''
        Imprime un pequeño mensaje de error.
        '''
        
        print("\033[1;31m   ,")
        print("  / \  ")
        print(" / ! \ ERROR:")
        print("/_____\  \033[00m\n")


def es_fasta (archivo):

	'''
	Comprueba que el archivo tiene formato fasta. En caso de no serlo, 
        mostrará el mensaje de ayuda y finalizará la ejecución del programa.
	'''

	with open(archivo, "r") as input_handle:

		 if not input_handle.read().startswith(">"):
                        error()
                        print(archivo+" no tiene formato fasta.\n")
                        ayudacorto()
                        sys.exit()


def es_prosite (archivo):

        """
        Comprueba que un archivo tiene formato Prosite. En caso de no serlo, 
        mostrará el mensaje de ayuda y finalizará la ejecución del programa.
        """

        with open(archivo, "r") as input_handle:

                if not input_handle.read().startswith("CC"):
                        error()
                        print(archivo+" no tiene formato Prosite.\n")
                        ayudacorto()
                        sys.exit()


def lista_gbk(listanombres,dirpath):

        """
        Crea una lista de archivos con formato GenBank.
        Si no tiene formato GenBank, no se analiza ese fichero.
        """

        gbk_list = []
        for nombre in listanombres:
                GBK = str(dirpath)+nombre
                gbk = open(str(dirpath)+nombre,'r').read()

                if gbk.startswith("LOCUS"): #Comprobamos que tenga formato gbk 
                        gbk_list.append(nombre)
                else: 
                        print("El archivo "+GBK+ " no tiene formato Genbank,"
                                +" por lo que no será analizado.")
        return gbk_list


def control_porcentaje(numero):

        """
        Comprueba si el argumento (coverage o identity) es un número 
        y si está entre 0 y 100. En caso de no serlo, mostrará el mensaje 
        de ayuda y finalizará la ejecución del programa.
        """

        #Creamos una variable para meter el número de forma provisional
        #sin ser decimal o negativo
        nprovisional = numero.replace("-", "").replace(".", "", 1)
        if not str(nprovisional).isdigit():
                error()
                print(numero+' no es un valor numérico.')
                ayudacorto()
                sys.exit()
        
        else:
                numero = float(numero)

                #Comprobamos que el número esté entre 0 y 100 
                #(ya que debe de ser un %) 
                if numero > 100:
                        numero = str(numero)
                        error()
                        print("El número introducido ("+numero+") es > 100.\n")
                        ayudacorto()
                        sys.exit()
                elif numero < 0:
                        numero = str(numero)
                        error()
                        print("El número introducido ("+numero+") es < 0.\n")
                        ayudacorto()
                        sys.exit()


def numero_argumentos():

        '''
        Detecta cuándo el número de argumentos no es correcto y muestra 
        el mensaje de ayudacorto. Además, finaliza el proceso del script.
        '''

        if len(sys.argv) != 6:
                error()
                print ("El número de argumentos es incorrecto.")
                print ("Argumentos a introducir: 5")
                ayudacorto()
                sys.exit()


def existe_archivo(archivo):

        '''
        Comprueba si el archivo existe. Si no, manda el mensaje de ayudacorto
        y sale del programa.
        '''

        if not os.path.isfile(archivo):
                error()
                print (archivo + " no existe.")
                ayudacorto()
                sys.exit()


def existe_dir(directorio):

        '''
        Comprueba si el directorio existe. Si no, manda el mensaje de ayudacorto
        y sale del programa.
        '''

        if not os.path.isdir(directorio):
                error()
                print(directorio + " no existe.")
                ayudacorto()
                sys.exit()


##########Código###########

def main():
        '''
        Función principal que ejecuta el código del script main.py

        DESCRIPCIÓN: Este script realiza las siguientes tareas:
                - BlastP de cada query aportada por el usuario 
                  frente a proteínas obtenidas de Genbank.
                - Alineamiento de los hits procedentes del BlastP 
                  (usando Muscle).
                - Árbol filogenético por el metodo Neighbor-Joining 
                  para cada proteína query (usando Muscle).
                - Búsqueda de dominios protéicos presentes en database 
                  de Prosite indicando sus características principales 
                  (name, accession, description, domain...).

        USO: ./main.py [query] [dir_subject/] [identity cut-off] 
                [coverage cut-off] [prosite database]
        '''        
        
        print(" ")

        #Antes de ejecutar el código principal, realizaremos los CONTROLES.

        #Si el primer argumento es -h o --help, se mostrará en pantalla el
        #mensaje de ayuda largo.
        ayuda()

        #Comprobamos que el número de argumentos sea correcto.
        numero_argumentos()

        #Definimos las variables de los argumentos.
        query = sys.argv[1]
        dir_subject = sys.argv[2]
        identity = sys.argv[3]
        coverage = sys.argv[4]
        prosite_db = sys.argv[5]

        #Comprobamos que el coverage y query sean numeros entre 0 y 100.
        control_porcentaje(coverage)  
        control_porcentaje(identity)

        #Comprobamos que el archivo query y el database existen.
        existe_archivo(query)
        existe_archivo(prosite_db)

        #Comprobamos que el directorio de los subject (gbk) existe.
        existe_dir(dir_subject)
        
        #Comprobamos que el archivo query tiene formato fasta 
        #y que prosite_database tiene formato prosite
        es_fasta(query)
        es_prosite(prosite_db)


        print ("\n\033[;36m"+"Por favor, tenga paciencia, el programa"
                +" se está ejecutando."+"\033[00m")


        #Creamos los directorios para almacenar datos, resultados...

        #Creamos un directorio principal con el nombre analisis_ junto con la
        #fecha (día-mes-año) y hora:minuto:segundo de creación del archivo. 
        now=datetime.now()
        fecha = (format(now.day)+"-"+format(now.month)+"-"+format(now.year)
                +"@"+format(now.hour)+":"+format(now.minute)+":"
                +format(now.second))
        os.makedirs("analisis_"+fecha+"/results")
        os.mkdir("analisis_"+fecha+"/data")
        os.mkdir("analisis_"+fecha+"/temporary_file")


        #Obtenemos una lista con los nombres de todos los subjects GBK para
        #hacer un único archivo multifasta. La función hacer_multifasta del 
        #módulo blast además comprueba que los archivos tengan forma gbk
        for (dirpath, dirnames, filenames) in os.walk(dir_subject):
                gbk_lista = lista_gbk(filenames,dirpath)

                for GBK in gbk_lista:
                        GBK = str(dirpath)+GBK
                        blast.hacer_multifasta(GBK, "subject_fasta")

        #Separamos cada query del archivo query para analizarlas por separado
        #y hacemos el análisis blast.
        with open(query, "r") as input_handle:

                for record in SeqIO.parse(input_handle, "fasta"):
                        sequence = str('>' + record.id + "\n" + record.seq)
                        archivo_temporal =  open("temporal0", "w")
                        archivo_temporal.write(sequence)
                        archivo_temporal.close()
                        blast.hacer_blast('temporal0', 'subject_fasta', record.id)
                        os.remove('temporal0')

                        #Añadimos la query a los resultados de BLAST filtrados
                        resultados_blast = open("resultados_blast_"+record.id,'a')
                        resultados_blast.write(str("\n"+sequence))
                        resultados_blast.close()

                                        
                        #Realizamos el árbol filogenético usando MUSCLE y el 
                        #alineamiento. En caso de que no obtengamos resultados
                        #en el blast tras filtrar (coverage e identity), no 
                        #realizaremos el alineamiento ni árbol de muscle. Para
                        #ello usamos try-except.          
                        try:
                                muscle.muscle("resultados_blast_"+record.id, 
                                        'arbol_muscle_'+record.id, record.id)
                        except:
                                pass
                                        

                        #Realizamos el análisis de PROSITE.
                        
                        #Obtenemos los subject completos (con query) para
                        #realizar el análisis.
                        subject_fasta = open("subject_fasta", 'a')
                        subject_fasta.write(str("\n"+sequence))
                        subject_fasta.close()  
                        blast.hacer_subject_prosite('subject_fasta',
                                'resultados_blast_'+record.id, record.id)

                        #Ejecutamos la función de prosite utilizando como
                        #argumentos el archivo fasta creado en el anterior
                        #paso, la database de prosite y almacenamos el 
                        #resultado en otro archivo (resultados_prosite_ID...) 
                        prosite_mod.prosite_f('subject_prosite_'+record.id, 
                                prosite_db, 'resultados_prosite_'+record.id)
                        
                        
                        #Creamos los directorios correspondientes para los 
                        #resultados en función de las querys.                        
                        os.mkdir("analisis_"+fecha+"/results" + "/"+record.id)
                        os.mkdir("analisis_"+fecha+"/results" + "/"+record.id
                                +"/blast")
                        os.mkdir("analisis_"+fecha+"/results" + "/"+record.id
                                +"/muscle")
                        os.mkdir("analisis_"+fecha+"/results" + "/"+record.id
                                +"/prosite")

                        #Copiamos los archivos que dependan de las querys a las
                        #rutas correspondientes.
                        shutil.move("resultados_blast_"+record.id, "analisis_"
                                +fecha+"/results"+"/"+record.id
                                +"/blast/resultados_blast_"+record.id)
                        shutil.move("provisional_"+record.id, "analisis_"
                                +fecha+"/temporary_file/provisional_"
                                +record.id)
                        shutil.move("alineamiento_"+record.id, "analisis_"
                                +fecha+"/results"+"/"+record.id
                                +"/muscle/alineamiento_"+record.id)
                        shutil.move("arbol_muscle_"+record.id, "analisis_"
                                +fecha+"/results"+"/"+record.id
                                +"/muscle/arbol_muscle_"+record.id)
                        shutil.move("resultados_prosite_"+record.id, "analisis_"
                                +fecha+"/results"+"/"+record.id
                                +"/prosite/resultados_prosite_"+record.id)
                        shutil.move("subject_prosite_"+record.id, "analisis_"
                                +fecha+"/data/subject_prosite_"+record.id)


        #Copiamos los archivos y dirsubject las rutas correspondientes para
        #que los ficheros estén ordenados.
        shutil.copytree(dir_subject, "analisis_"+fecha+"/data/"+dir_subject)
        shutil.copy(prosite_db, "analisis_"+fecha+"/data/"+prosite_db)
        shutil.copy(query, "analisis_"+fecha+"/data/"+query)
        shutil.move("subject_fasta", "analisis_"+fecha+"/data/subject_fasta")

        #Sé que esta última parte de los directorios es un poco liosa pero llevo dos semanas 
        #durmiendo muy poco por las recus, lo siento :(


        print("\033[1;36m"+"\n¡El proceso ha finalizado!\n")
        print("\033[0;36m"+"\nLos resultados y datos se han almacenado en "
                +"la carpeta: analisis_"+fecha+"\n")


#Ejecutamos el código con la función principal
main()