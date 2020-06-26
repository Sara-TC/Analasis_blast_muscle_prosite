#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########Módulos##########

import sys
import re
from Bio.ExPASy import Prosite,Prodoc
from Bio import Seq
from Bio import SeqIO

#########Funciones#########

def prosite_f (query_file, domains_file, output):

    '''
    Encuentra los dominios de una proteína (presente en el archivo query_file)
    comparando dominios de una database de prosite (domains_file) en formato 
    prosite. La información del análisis se almacena en el archivo output.
	'''

    with open(query_file,"r") as handle:

            for record_fasta in SeqIO.parse(handle, "fasta"):

                out = open(output, 'a')
                with open(domains_file) as handle2:

                    #Parsear los Records y recorrerlos
                    for record3 in Prosite.parse(handle2):

                        #Obtenemos un string con la secuencia del query
                        query = str(record_fasta.seq)

                        #Obtenemos cada dominio de prosite traducido 
                        #(a lenguaje python)
                        dom = str(record3.pattern.replace("}","]"
                            ).replace("{","[^").replace(".", "").replace("x", "."
                            ).replace(")", "}").replace("(", "{").replace("-", ""))

                        #Comprobamos si el dominio aparece en la secuencia query
                        if re.search(dom, query) and dom != "":

                            #Escribir la información en el fichero output
                            out.write("-Blast hit: "+record_fasta.id+"\n")
                            out.write("Name: "+record3.name+"\n")
                            out.write("Accession: "+record3.accession+"\n")
                            out.write("Description: "+record3.description+"\n")
                            out.write("Pattern: "+dom+"\n")
                            out.write("Domain: "+str(re.findall(dom, query)[0])
                                    +"\n\n")

                #Cerrar el fichero de resultados
                out.close()

####################