#!/usr/bin/python
# -*- coding: utf-8 -*-

#transforme le fichier fichier_dep en fichier_arr avec les couleurs et les tailles de fichier_info (ne fait pas de calcul de taille)
#le format des fichiers dep et arr, sont des sorties de cellDesigner
#fichier info a une ligne par reaction a modifier : nom (tabulation) taille (tabulation) code couleur (attention prendre en compte le premier champ (qu'est-ce que c'est ? dans le doute mettre ff puis le code hexa de la couleur voulue)
#import pickle


import pandas

def efm_sbml(fichier_dep,fichier_info,fichier_arr):
	import re
	dico=parserInfo(fichier_info)
	
	file_dep=open(fichier_dep,'r')
	lines=file_dep.readlines()
	file_dep.close()
	
	file_arr=open(fichier_arr,'w')
	inReac=False
	for line in lines:
		if inReac : #on est dans une reaction a modifier
			if line[0:19]=="<celldesigner:line ":
				#on modifie
				regex='width="(.+?)"'
				line=re.sub(regex,'width="'+dico[name]["width"]+'"',line)
				regex='color="(.+?)"'
				line=re.sub(regex,'color="'+dico[name]["color"]+'"',line)
			elif line=="</annotation>\n":
				inReac=False
			#on ecrit
			file_arr.write(line)
		else :
			#on recopie la ligne
			file_arr.write(line)
			if line[0:9]=="<reaction" : #on est dans une reaction, a modifier ?
				#recuperer name
				regex='name="(.+?)"'
				search=re.search(regex,line)
				name=search.group(1)
				#voir si la reaction est a modifier
				if name in dico:
					inReac=True
	file_arr.close()
	
def parserInfo(lines):
	dico={}
	#file_=open(fichier,'r')
	#lines=file_.readlines()
	#file_.close()
	
	for line in lines :
		dico2={}
		tab=line.split()
		dico2["width"]="4" #tab[1]
		dico2["color"]=tab[2] #"ffff0000"
		dico[tab[0]]=dico2
	return dico

def pkl_read(xmlfile, pklfile):

	#with open(pklfile, 'rb') as f:
		#df = pickle.load(f) # df = dataframe Pandas
	df = pandas.read_csv(pklfile, index_col=0)
	for i, efm in df.iterrows(): # liste de efms
		lines = []
		supp = efm.loc[efm != 0] # support
		for index, value in supp.iteritems():
			color = 'ff0000ff' if value >= 0 else 'ff0000ff'#'ff00ff00'
			line = str(index) + ' 5.0 ' + color + '\n'
			lines.append(line)
		efm_sbml(xmlfile, lines, xmlfile[:-3]+f'_efm_{i}'+'.xml')
	
if __name__ == "__main__":
	import sys
	pkl_read(sys.argv[1], sys.argv[2])
