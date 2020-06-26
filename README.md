# Analasis_blast_muscle_prosite

#### DESCRIPCIÓN 

Este script realiza las siguientes tareas:

        - BlastP de cada query aportada por el usuario frente a proteínas obtenidas de Genbank.
        - Alineamiento de los hits procedentes del BlastP (usando Muscle).
        - Árbol filogenético por el metodo Neighbor-Joining para cada proteína query (usando Muscle).
        - Búsqueda de dominios protéicos presentes en database de Prosite indicando sus
          características principales (name, accession, description, domain...).
          
          
#### INSTALACIÓN

Muscle y Biopython.


#### USO 

    ./main.py [query] [dir_subject/] [identity cut-off] [coverage cut-off] [prosite database]
    
   *[query]*: archivo a analizar con formato fasta.
  
   *[dir_subject/]*: directorio donde se encuentran los subjects con formato Genbank.
   
   *[identity cut-off]*: número entre 0 y 100.
   
   *[coverage cut-off]*: número entre 0 y 100.
   
   *[prosite database]*: base de datos de Prosite con formato Prosite.
   
   
    
    ./main.py [-h]
    
   o
   
    ./main.py [--help]
    
   Mostrará un mensaje de ayuda.    
