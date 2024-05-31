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
import checkResults


 #------------------EFM SOLVING----------------
def efm_solving(model,efm,dico_param):
    
    #CREATION OF THE SUBMODEL ASSOCIATED TO THE EFM (CONTAINS ONLY THE METABOLITES AND REACTIONS IN THE EFM)
    submodel = submodelCreation.submodel(model,efm)

    #W CONTAINING ONLY THE REACTIONS OF THE EFM
    W_efm = dico_param["W"].iloc[:,function.find_efm(efm)]
    
    #DICTIONARY CONTAINING THE INDICES OF CERTAIN METABOLITES IN THE MODEL 
    dico_metabolites_indexes = submodelCreation.metabolites_indexes(submodel)

    #EFM'S RESOLUTION WITH CVXYPY
    result = cvxEfm.cvx_problem(submodel,efm,dico_param,model,W_efm,dico_metabolites_indexes) 

    return result


 #------------------STORAGE OF RESULTS IN AN EXCEL FILE-----------------
def result_file_redaction(model,efm,result_cvx,efm_number,constraints,dG,ratio,redox,ph):

    #EXTRACTION OF CVX RESULTS
    flux  = []
    for value in efm :
        flux.append(float(result_cvx["alph"])*float(value)*3600) # flux en h-1
    
    
    enzyme = np.transpose(result_cvx["Enz"])
    concentration = np.transpose(result_cvx["X"])
    metabolite = [model.metabolites]
    reactions = [model.react_name]
    
    #WRITING THE RESULTS FILE
    workbook = xlwt.Workbook()
    function.excel_file(workbook,"Metabolites",["Metabolite", "Concentration"],metabolite,concentration)
    function.excel_file(workbook,"Enzymes",["Reaction","Flux (mmol.gDW-1.h-1)","Enzyme"],reactions,flux,enzyme)
    function.excel_file(workbook,"Constraints",["Reaction","contraintes"],reactions,constraints)
    function.excel_file(workbook,"DeltaG",["Reaction","dG"],reactions,dG)
    function.excel_file(workbook,"Ratio",["Redox","Value"],ratio)
    function.excel_file(workbook,"Redox",["Energetic potentials","Value"],redox)
    function.excel_file(workbook,"pH",["pH","Value"],ph)
    workbook.save('results/results_efm_'+str(efm_number)+'.xls')

    return 0


#------------------APPLICATION OF "efm_solving" AND "result_file" ON THE EFM(s) OF INTEREST-----------------
def efm_definition(efm_number,model,dico_param,outputfile):

    #INDICES OF THE EFM(s) TO SOLVE AMONG THOSE HAVING A BIOMASS DIFFERENT FROM 0
    k_opt = [] 
    efm_name = []
    
    for i in range(len(dico_param["lobj"])):
        if efm_number == "all":
            efm_name.append(model.efms.axes[0][dico_param["lobj"][i]]) 
            k_opt.append(i) 
        else: #If only a particular efm to solve, storage of its number and indexe 
            if model.efms.axes[0][dico_param["lobj"][i]] == int(efm_number) :
                efm_name.append(model.efms.axes[0][dico_param["lobj"][i]])
                k_opt.append(i)

    #RESOLUTION OF THE ENZYMATIC PROFILE OF THESE EFM(s)
    indice = 0
    dico_opt = {}
    dico_enz = {}

    for efm_indexe in k_opt: 
        efm = dico_param["EM"].iloc[dico_param["lobj"][efm_indexe],:] #storage of the efm to solve 
        for i in efm :
            i = float(i) 
        print("efm :",efm_name[indice]) #display of the efm solved 
        try:
            ko2 = function.indice_find(model.reactions,'ex-o2',"name","")
            cvx_result = efm_solving(model,efm,dico_param) #efms solving with CVXPY

            dico_opt[efm_name[indice]]= cvx_result['opt'] 
            dico_enz[efm_name[indice]]= np.sum(cvx_result['Enz'])
        except cvxpy.error.SolverError as e:
            dico_opt[efm_name[indice]]= -1
            dico_enz[efm_name[indice]]= -1
        # excel output    
        constraints,dG,ratio,redox,ph=checkResults.check(efm,dico_param,model,cvx_result)
        result_file_redaction(model,efm,cvx_result,efm_name[indice],constraints,dG, ratio,redox,ph) #redaction of the results file 
        
        indice += 1                
        

    dico_opt= dict(sorted(dico_opt.items(), key=operator.itemgetter(1),reverse=True))
    print(dico_opt)
    with open(outputfile,'w') as f:
        for key, value in dico_opt.items(): 
            f.write('%s : %s : %s\n' % (key, value, dico_enz[int(key)]))

   

    

    return 0





