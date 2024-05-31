import numpy as np
import cvxpy as cp
from math import *
import basicFunctions as function

#--------------RESOLUTION OF THE ENZYMATIC PROFILE-------------- 
def cvx_problem(submodel_dico,efm,dico_params,model,W_efm,dico_ind):
    
    ##DEFINITIONS OF THE PROBLEM VARIABLES
    x = cp.Variable(len(submodel_dico["submatrix_metab"]),pos=True) #CVX variable containing the concentration of all the efm metabolites 
    alph = cp.Variable(pos=True) #CVX variable to resolve the flux rate of each reaction
    E = cp.Variable((len(function.find_efm(efm)),1),pos=True) #CVX variable containing the enzymatic concentration of the flows 
   
    #print(submodel_dico["submatrix_metab"])
    #TURN C AND e INTO VECTORS TO MULTIPLY THEM IN THE OBJECTIVE 
    c_vector = function.create_vector(dico_params["C"].loc[1,:]) #Turns the column of C containing the numerical values into a vector
    e_vector = function.create_vector(efm) #Turns e into a vector
    

    #DEFINITIONS OF THE PROBLEM OBJECTIVE
    objective = cp.Maximize(alph*(c_vector@e_vector))


    #DEFINITIONS OF THE PROBLEM CONSTRAINTS
    constraints = []
    
    i = 0 #Index of the enzymes of E
    j = 0 #Index of efm reactions

    while j < len(efm) : #For each reaction of the efm 

        

        if efm[j] != 0 : #if the reaction value is different from 0

            lsub_npArray = function.find(submodel_dico["submatrix"][:,j]<0) #Index of reaction substrates in the submodel matrix 
            lprod_npArray = function.find(submodel_dico["submatrix"][:,j]>0) #Index of reaction products in the submodel matrix 
            lsubstrat_npArray = function.find(model.matrice[:,j]<0) #Index of reaction substrates in the model matrix 
            lproduit_npArray = function.find(model.matrice[:,j]>0) #Index of reaction products in the model matrix 

            #If it is a diffusion reaction, recuperation of the associated diffusion coefficient
            idxDiff = function.dData(model,model.diffData,j) 

            #Applications of diffusion constraints
            if idxDiff != []:

                if model.react_name[j] == 'ex-o2':
                    constraints.append(x[dico_ind["k_O2"]]==0.21)
                    constraints.append(x[dico_ind["k_OXY_ext"]]==0.21)
                if model.react_name[j] == 'ex-co2':    
                    constraints.append(x[dico_ind["k_CO2"]]==0.34)
                    constraints.append(x[dico_ind["k_CO2_ext"]]==0.34)
                if model.react_name[j] == 'ex-etoh':
                    constraints.append(x[dico_ind["k_ETOH"]]==109)
                    constraints.append(x[dico_ind["k_ETOH_ext"]]==109)

                constraints.append(E[i]<=1e-12) 
                i+=1
                j+=1
            

            #If it is not an diffusion reaction 
            else : 
                #constraints.append(E[i]>=1e-10)
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
                    #constraints.append(E[i]>=1e-4)
                    sH = 4
                    zetak=0.86
                    zetaP=0.14
                    gammak=1
                    gammaP=0.419
                   
                else : 
                    sH = 0
                
                #If it is the biomass reaction
                if j == function.find(dico_params["C"].iloc[1,:])[0] : 
                    T1=function.produitBiom(x,lsub_npArray,lsubstrat_npArray,model.km,j)
                    T2=function.produitBiom(x,lprod_npArray,lproduit_npArray,model.km,j)
                    #constraints.append(((alph*efm[j]))/(model.kcatp.iloc[j,1]*E[i]) <=1) 
                    i+=1
                    j+=1

                #For the other reactions
                else :
                    T1=function.produit(x,lsub_npArray,lsubstrat_npArray,model.km,submodel_dico,j)
                    T2=function.produit(x,lprod_npArray,lproduit_npArray,model.km,submodel_dico,j)
                    Hint=model.constantes['H_int']
                    a=model.constantes['a']
                    b=model.constantes['b']
                    F=model.constantes['F']
                    R=model.constantes['R']
                    T=model.constantes['T']

                    if float(efm[j]) > 0 :
                        
                        if sH != 0:
                            W1=(x[dico_ind["k_Hext"]]/Hint)**(-sH*gammak*a*F/(log(10)*R*T)) * exp(-sH*gammak*b*F/(R*T))
                            W2=(x[dico_ind["k_Hext"]]/Hint)**((1-gammak)*sH*a*F/(log(10)*R*T)) * exp((1-gammak)*sH*b*F/(R*T))
                            constraints.append((((1+T1+T2)*alph*efm[j]) + (model.kcatm.iloc[j,1]*T2*E[i]*W2))/(model.kcatp.iloc[j,1]*T1*E[i]*W1)<=1)
                            
                        else :
                            constraints.append((((1+T1+T2)*alph*efm[j]) + (model.kcatm.iloc[j,1]*T2*E[i]))/(model.kcatp.iloc[j,1]*T1*E[i])<=1)
                        i+=1
                        j+=1

                    elif float(efm[j]) < 0:
                        
                        if sH != 0 :
                            W1=(x[dico_ind["k_Hext"]]/Hint)**(-sH*gammak*a*F/(log(10)*R*T)) * exp(-sH*gammak*b*F/(R*T))
                            W2=(x[dico_ind["k_Hext"]]/Hint)**((1-gammak)*sH*a*F/(log(10)*R*T)) * exp((1-gammak)*sH*b*F/(R*T))
                            constraints.append((((1+T1+T2)*alph*abs(efm[j])) + (model.kcatp.iloc[j,1]*T1*E[i]*W1))/(model.kcatm.iloc[j,1]*T2*E[i]*W2)<=1)
                        else :
                            constraints.append((((1+T1+T2)*alph*abs(efm[j])) + (model.kcatp.iloc[j,1]*T1*E[i]))/(model.kcatm.iloc[j,1]*T2*E[i])<=1)
                        i+=1
                        j+=1  
        else: 
            j+=1  


    #CONSTRAINTS ON METABOLITES 
    #CONSTRAINTS ON INTERNAL METABOLITES
    idx=len(submodel_dico["submatrix_metab_sorted"][0])+1
    sumX=0
    for k in submodel_dico["submatrix_metabInt_idx"]:
        if submodel_dico["submatrix_metab"][k].name != 'H':
            constraints.append(x[k] >= 1e-4) # contrainte tEFMA cmin>= 10-7 M
        sumX=sumX+x[k]
    constraints.append(sumX <= 300) # ref Bennet 




    #TURN W AND E INTO VECTORS TO MULTIPLY THEM IN A CONSTRAINT
    w_vector = function.create_vector( W_efm.iloc[1,:])
    E_vector = function.create_vector(E)

    constraints.append(w_vector@E_vector<= model.constantes['Etot'])

   
    constraints.append((x[dico_ind["k_H"]] == model.constantes['H_int']))
    constraints.append(x[dico_ind["k_Hext"]]*1e-3 <= exp(-log(10)*model.constantes['pHMin']))
    constraints.append(exp(-log(10)*model.constantes['pHMax']) <= x[dico_ind["k_Hext"]]*1e-3 )
    


   #CONSTRAINTS ON EXTERNAL METABOLITES
   
   #STORAGE OF THE SUBMODEL EXTERNAL METABOLITE INDICES IN THE MODEL
    listofmetabext = []
    ind = 0 
    for met in model.ext_met :
        for metSubmat in submodel_dico["submatrix_metab_sorted"][1]:
            if met == metSubmat:
                listofmetabext.append(ind)
               
        ind += 1

    #CREATION OF A LIST l CONTAINING THE VALUES OF EXTERNAL METABOLITES SORTED ACCORDING TO THE POSITION OF THE EXTERNAL METABOLITES IN THE MODEL
    listmetextsorted = []
    for i in model.xext.iloc[:,0]:
        listmetextsorted.append(0)

    metInd = 0
    for i in model.xext.iloc[:,0]:
        ind = 0
        for j in model.ext_met:
            if i == j.name :
                listmetextsorted[ind] = model.xext.iloc[metInd,1]
            ind += 1
        metInd += 1


    ind = 0   
    for k in range(len(submodel_dico["submatrix_metab_sorted"][1])):
        constraints.append(x[k] <= listmetextsorted[listofmetabext[ind]])
                
        ind += 1
        
   
    #CVX PROBLEM DEFINITION AND RESOLUTION
    problem = cp.Problem(objective, constraints)
    problem.solve(gp=True)

    #PROBLEM SATUTS
    print("=================")
    print(problem.status)
    print("=================")

    opt = problem.value * 3600
    cvxStatus = problem.status

    
    
    #STORAGE OF THE CVX VARIABLES VALUES
    lk = function.find_efm(efm)
    Enz = function.remplissage_val(efm,lk,E)
    

    X = function.remplissage_val(model.metabolites,submodel_dico["submatrix_metab_idx"],x)

    dico_result = {"opt":opt,"X":X,"Enz":Enz,"alph":alph.value,"cvxStatus":cvxStatus}
    
    return dico_result


