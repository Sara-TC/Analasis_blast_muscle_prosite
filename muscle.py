#!/usr/bin/env python
# -*- coding: UTF-8 -*-

##########Módulos##########

import os

#########Funciones#########

def muscle(input, output, ID):

    '''
    Alinea un archivo fasta (input) usando el comando de muscle 
    (perteneciente a bash) y realiza un árbol filogenético por el 
    método Neigbour Joining, el cual se almacena en un archivo (output).
    '''
    
    #usamos 2>/dev/null para que no nos aparezcan en pantalla el output de muscle
    x = "muscle -in "+input+" -out alineamiento_"+ID+" 2>/dev/null" 
    os.system(x)
    
    y = "muscle -maketree -in "+'alineamiento_'+ID+" -out "+output+" -cluster neighborjoining" + " 2>/dev/null"
    os.system(y)

##########################