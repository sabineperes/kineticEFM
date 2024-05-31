import numpy as np
import cvxpy
from math import *
import basicFunctions as function
import cvxEfm 
import submodelCreation
import xlwt
import asp_results 
import pprint
import operator

def metab_indexes(model):

    k_H=function.indice_find(model.metabolites,'H',"name","")
    k_Hext=function.indice_find(model.metabolites,'H-ext',"name","")
    k_ATP=function.indice_find(model.metabolites,'ATP',"name","")
    k_ADP=function.indice_find(model.metabolites,'ADP',"name","")
    k_NAD=function.indice_find(model.metabolites,'NAD',"name","")
    k_NADH=function.indice_find(model.metabolites,'NADH',"name","")
    k_NADP=function.indice_find(model.metabolites,'NADP',"name","")
    k_NADPH=function.indice_find(model.metabolites,'NADPH',"name","")
    k_Q8=function.indice_find(model.metabolites,'Q8',"name","")
    k_Q8H2=function.indice_find(model.metabolites,'Q8H2',"name","")
    k_ACOA=function.indice_find(model.metabolites,'ACETYL-CoA',"name","")
    k_COA=function.indice_find(model.metabolites,'CoASH',"name","")
    k_SUCC_COA=function.indice_find(model.metabolites,'SUCC-CoA',"name","")
    k_O2=function.indice_find(model.metabolites,'O2',"name","")
    k_OXY_ext = function.indice_find(model.metabolites,'OXY-ext',"name","")
    k_CO2 = function.indice_find(model.metabolites,'CO2',"name","")
    k_CO2_ext = function.indice_find(model.metabolites,'CO2-ext',"name","")
    k_ETOH = function.indice_find(model.metabolites,'ETOH',"name","")
    k_ETOH_ext = function.indice_find(model.metabolites,'ETOH-ext',"name","")

    dico = {"k_H":k_H, "k_Hext":k_Hext, "k_ATP":k_ATP, "k_ADP":k_ADP, "k_NAD":k_NAD, "k_NADH":k_NADH, "k_NADP":k_NADP, "k_NADPH":k_NADPH,"k_Q8":k_Q8, "k_Q8H2":k_Q8H2, "k_ACOA":k_ACOA, "k_COA":k_COA, "k_SUCC_COA":k_SUCC_COA, "k_O2":k_O2, "k_OXY_ext":k_OXY_ext, "k_CO2":k_CO2, "k_CO2_ext":k_CO2_ext, "k_ETOH":k_ETOH, "k_ETOH_ext":k_ETOH_ext}
   
    return dico

def check(efm,dico_params,model,dico_result):


    #CREATION OF THE SUBMODEL ASSOCIATED TO THE EFM (CONTAINS ONLY THE METABOLITES AND REACTIONS IN THE EFM)
    submodel = submodelCreation.submodel(model,efm)
    #W CONTAINING ONLY THE REACTIONS OF THE EFM
    W_efm = dico_params["W"].iloc[:,function.find_efm(efm)]
    #DICTIONARY CONTAINING THE INDICES OF CERTAIN METABOLITES IN THE MODEL 
    dico_ind = metab_indexes(model)
    
    ##DEFINITIONS OF THE PROBLEM VARIABLES
   
    E = np.transpose(dico_result["Enz"])
    x = np.transpose((dico_result["X"]))
    alph = dico_result["alph"]

    constraints = []
    dG = []


    pHint = -log10(x[dico_ind["k_H"]]*1e-3)
    pHext = -log10(x[dico_ind["k_Hext"]]*1e-3)
    deltapH = pHint - pHext
    deltaPhi = model.constantes['a']*deltapH + model.constantes['b']
    pmf= -deltaPhi + (log(10) * model.constantes['R'] *model.constantes['T'] * deltapH)/model.constantes['F']
    print('PMF')
    print(pmf)
    
    if x[dico_ind["k_O2"]]==0:
        ratio=nan
    else:
        ratio=float(x[dico_ind["k_H"]]*2/x[dico_ind["k_O2"]])

    dico_pH={"pHint": pHint, "pHext": pHext}
    dico_redox={"PMF": pmf, "nadh/nad": float(x[dico_ind["k_NADH"]]/x[dico_ind["k_NAD"]]),"nadph/nadp": float(x[dico_ind["k_NADPH"]]/x[dico_ind["k_NADP"]]), "atp/adp": float(x[dico_ind["k_ATP"]]/x[dico_ind["k_ADP"]]), "H/O": ratio }
    dico_ratio={"pHint": pHint, "pHext": pHext, "deltapH": deltapH, "deltaPsi": deltaPhi, "PMF": pmf, "nad/nadh": float(x[dico_ind["k_NAD"]]/x[dico_ind["k_NADH"]]), \
    "nadp/nadph": float(x[dico_ind["k_NADP"]]/x[dico_ind["k_NADPH"]]), "atp/adp": float(x[dico_ind["k_ATP"]]/x[dico_ind["k_ADP"]]), "H/O2": ratio}
    
    i = 0 #Index of the enzymes of E
    j = 0 #Index of efm reactions
    

    while j < len(efm) : #For each reaction of the efm 
        if efm[j] != 0 : #if the reaction value is different from 0

            
            lsubstrat_npArray = function.find(model.matrice[:,j]<0) #Index of reaction substrates in the model matrix 
            lproduit_npArray = function.find(model.matrice[:,j]>0) #Index of reaction products in the model matrix 
            
            #If it is a diffusion reaction, recuperation of the associated diffusion coefficient
            idxDiff = function.dData(model,model.diffData,j) 

            #Applications of diffusion constraints
            if idxDiff != []:
                if model.react_name[j] == 'ex-o2':
                    constraints.append(1)#x[dico_ind["k_O2"]])
                if model.react_name[j] == 'ex-co2':    
                    constraints.append(1)#x[dico_ind["k_CO2"]])
                if model.react_name[j] == 'ex-etoh':
                    constraints.append(1)#x[dico_ind["k_ETOH"]])
                dG.append(-1)
                 
                i+=1
                j+=1
            

            #If it is not an diffusion reaction 
            else : 

                if model.react_name[j] == 'nadh-dh':
                    sH = 4
                    zetak=0.8
                    zetaN=0.17
                    zetaQ=0.009
                    gammak=0.049
                    gammaN=1
                    gammaQ=0.997
                    
                elif model.react_name[j] == 'cytbd':
                    sH = 2
                    zetak=0.477
                    zetaQ=0.02
                    zetaO=0.5
                    gammak=0.037
                    gammaQ=0.408
                    gammaO=0.39
                   
                elif model.react_name[j] == 'atp-synth' : 
                    sH = 4
                    zetak=0.86
                    zetaP=0.14
                    gammak=1
                    gammaP=0.419
                   
                else : 
                    sH = 0
                

                #If it is the biomass reaction => nothing
                if j == function.find(dico_params["C"].iloc[1,:])[0] : 
                    dG.append(-1)
                    constraints.append(0)
                    #constraints.append(((alph*efm[j]))/(model.kcatp.iloc[j,1]*E[i]) <=1)
                    i+=1
                    j+=1

                #For the other reactions
                else :
                        
                    T1=function.produit2(x,lsubstrat_npArray,model.km,model,j)
                    T2=function.produit2(x,lproduit_npArray,model.km,model,j)

                    
                    TT1=function.produitdG(x,lsubstrat_npArray, model.matrice,j)
                    TT2=function.produitdG(x,lproduit_npArray,model.matrice,j)
                    

                    Hint=model.constantes['H_int']
                    a=model.constantes['a']
                    b=model.constantes['b']
                    F=model.constantes['F']
                    R=model.constantes['R']
                    T=model.constantes['T']

                    if float(efm[j]) > 0 :
                       
                        if sH != 0 :
                            W1=(x[dico_ind["k_Hext"]]/Hint)**(-sH*gammak*a*F/(log(10)*R*T)) * exp(-sH*gammak*b*F/(R*T))
                            W2=(x[dico_ind["k_Hext"]]/Hint)**((1-gammak)*sH*a*F/(log(10)*R*T)) * exp((1-gammak)*sH*b*F/(R*T))
                            constraints.append((((1+T1+T2)*alph*efm[j]) + (model.kcatm.iloc[j,1]*T2*E[i]*W2))/(model.kcatp.iloc[j,1]*T1*E[i]*W1))
                            dG.append(model.constantes['R']*model.constantes['T']*(log(TT2/TT1)-log(model.keq.iloc[j,1]))  +(model.constantes['F']*deltaPhi)- (sH*model.constantes['R']*model.constantes['T']*log(10)* deltapH))
                            
                        else :
                            constraints.append((((1+T1+T2)*alph*efm[j]) + (model.kcatm.iloc[j,1]*T2*E[i]))/(model.kcatp.iloc[j,1]*T1*E[i]))
                            dG.append(model.constantes['R']*model.constantes['T']*(log(TT2/TT1)-log(model.keq.iloc[j,1])))
                            
                           
                        i+=1
                        j+=1

                    elif float(efm[j]) < 0: #### verifier deltaG !!!!!!!!!!!!!!!!!
                        print('rev')
                        
                        if sH != 0 :

                            W1=(x[dico_ind["k_Hext"]]/Hint)**(-sH*gammak*a*F/(log(10)*R*T)) * exp(-sH*gammak*b*F/(R*T))
                            W2=(x[dico_ind["k_Hext"]]/Hint)**((1-gammak)*sH*a*F/(log(10)*R*T)) * exp((1-gammak)*sH*b*F/(R*T))
                            constraints.append((((1+T1+T2)*alph*abs(efm[j])) + (model.kcatp.iloc[j,1]*T1*E[i]*W1))/(model.kcatm.iloc[j,1]*T2*E[i]*W2))
                            dG.append(model.constantes['R']*model.constantes['T']*(log(TT1/TT2)+log(model.keq.iloc[j,1])) + (model.constantes['F']*deltaPhi)- (sH*model.constantes['R']*model.constantes['T']*log(10)* deltapH))
                           
                        else :
                            constraints.append((((1+T1+T2)*alph*abs(efm[j])) + (model.kcatp.iloc[j,1]*T1*E[i]))/(model.kcatm.iloc[j,1]*T2*E[i]))
                            dG.append(model.constantes['R']*model.constantes['T']*(log(TT1/TT2)+log(model.keq.iloc[j,1])))
                            
                        i+=1
                        j+=1  
        else: 
            j+=1  
            i+=1
            constraints.append(0)
            dG.append(nan)
    
    
    
    
    
    #dico_result = {"opt":opt,"X":X,"Enz":Enz,"alph":alph.value,"cvxStatus":cvxStatus}
    
    return constraints,dG,dico_ratio,dico_redox,dico_pH