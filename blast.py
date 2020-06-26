#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########Módulos##########

import sys
from Bio import Seq
from Bio import SeqIO
import re
import os

#########Variables#########

#Debemos escribir try except para que en caso de que se ejecute s
# main.py -h o --help no aparezca un error list index out of range
try: 
	query = sys.argv[1]
	dir_subject = sys.argv[2]
	identity = sys.argv[3]
	coverage = sys.argv[4]
except:
	pass

#########Funciones#########

def hacer_multifasta (subject, multifasta):

	"""
	Convierte un archivo con formato Genbank en formato multifasta.
	Además, añade un header con formato Locus_tag@Species_name.
	"""

	with open(subject, "r") as input_handle:

		for record in SeqIO.parse(input_handle, "genbank"):
			species = record.annotations["organism"].replace(" ","_")
						
			for feature in record.features:

				if feature.type == 'CDS':
					sys.stdout = open(multifasta, 'a')

					try:
						if feature.qualifiers['translation'][0]:
							print (">"+feature.qualifiers['locus_tag'][0]+"@"+species)
							print (feature.qualifiers['translation'][0])
					except:
						pass

		sys.stdout.close()
		sys.stdout = open("/dev/stdout",'w')


def hacer_blast (query, subject_fasta, ID):

	"""
	Realiza blast del query contra el subject aplicando los 
	criterios del identity y coverage cut-off.
	"""

	#Realiza el analisis de blast	
	analisis = "blastp -query " + query + " -subject " + subject_fasta + " -evalue '0.00001' -outfmt '6 qseqid sseqid qcovs pident evalue sseq' -out provisional_"+ID
	os.system(analisis)

	#Filtramos los resultados del blast	aplicando los criteros de coverage e indentity
	p = open('provisional_'+ID,"r")
	for line in p.readlines():
		fields = line.rstrip().split("\t")
		ide = fields[3]
		cov = fields[2]

		if ide >= identity and cov >= coverage:
			f=open('resultados_blast_'+ID,'a') #Los resultados del blast filtrados se almacenan en resultados_blast_<ID de cada query> 
			f.write('>'+fields[1]+'\n'+fields[5]+'\n')
			f.close()

	p.close()


def hacer_subject_prosite(subject_completos, subject_blast, ID):

	'''
	Genera unos subjects con secuencias completas para el análisis PROSITE.
	Para ello compara los nombres de los subjects tras el filtrado de blast
	(subject_blast) y almacena las secuencias en un nuevo fichero 
	(subject_prosite_ID). 
	'''

	sp=open("subject_prosite_"+ID, 'a')

	with open (subject_completos, 'r') as multifasta:

		for n in SeqIO.parse(multifasta, "fasta"):

			with open (subject_blast, 'r') as handle:

				for m in SeqIO.parse(handle, "fasta"):
					
					if str(m.id)==str(n.id):
						sp.write(">"+str(n.id)+"\n"+str(n.seq)+"\n")

	sp.close()

####################