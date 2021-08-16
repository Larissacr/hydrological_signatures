#------------Script for Hydrological Indices-------------#
# --Instituto de Pesquisas Hidraulicas - Hidrologia em Grande Escala (HGE)
# --Author: Larissa de Castro Ribeiro (larissa.ribeirocr@gmail.com)
# --Professor Advison: Dr. Rodrigo Cauduro Dias de Paiva
# --Date: 19/08/2020
# --Update in: 20/07/2021


# This script calculates de following hydrological indices:
    
#1 - Coefficient of baseflow recession (k)
#   Beck, H. E., van Dijk, A. I. J. M., Miralles, D. G., de Jeu, R. A. M., 
#   Sampurno Bruijnzeel, L. A., McVicar, T. R., & Schellekens, J. 2013. 
#   Global patterns in base flow index and recession based on streamflow 
#   observations from 3394 catchments. Water Resources Research, 49(12), 
#   7843â7863. doi:10.1002/2013wr013918.

#2 - Index of Aridity (IA)
#   Penman, A. D. M. 1953. Shear characteristics of a saturated silt, measured 
#   in triaxial compression. Geotechnique, 3(8), 312-328.

#3 - Hydrogram assymetry index (S)
#   Fleischmann, A. S., Paiva, R. C. D., Collischonn, W., Sorribas, M. V., & 
#   Pontes, P. R. M. (2016). On river-floodplain interaction and hydrograph
#   skewness. Water Resources Research, 52(10), 7615â7630. doi:10.1002/2016wr019233 


#4 - Baseflow Index (BFI)
#   Collischonn, Walter, and Fan, F.M., 2013, Defining parameters for Eckhardt’s 
#   digital baseflow filter: Hydrological 
#   Processes, v. 27, no. 18, p. 2,614–2,622, doi:10.1002/hyp.9391, 
#   http://dx.doi.org/10.1002/hyp.9391.

#5 – Streamflow Elasticity (STR_ela)
#   Sawicz, K., Wagener, T., Sivapalan, M., Troch, P. A., and Carrillo, G
#   Catchment classification: empirical analysis of hydrologic similarity based on
#   catchment function in the eastern USA. Hydrol. Earth Syst. Sci., 15, 2895–2911, 2011.
#   https://doi.org/10.5194/hess-15-2895-2011.

#6 - Maximum standardized reference flow (Q10_norm)
#   Sawicz, K. A., Kelleher, C., Wagener, T., Troch, P., Sivapalan, M., & Carrillo, G.
#   Characterizing hydrologic change through catchment classification.
#   Hydrology & Earth System Sciences Discussions, 10(5). 2014.
#   https://doi.org/10.5194/hess-18-273-2014.




########################### IMPORTING LIBRARIES ###########################

import pandas as pd
import datetime
import numpy as np
import winsound
from scipy.io import loadmat

#################################################################     
###                  SECTION 1 - IMPORTING DISCHARGE DATA                   ###
##################################################################
# The inputs for this script are: discharge, precipitation and/or evapotranspiration
# in DataFrama format. Lines represent different entries with a fixed time interval
# and the column are river segments.


###      Transforming input file of discharge from .mat to DataFrame         ###
##     Sub-catchments/River segments
mat1 = loadmat("D:/working_directory/discharge.mat")
mdata1 = mat1['...']
df_q = 0
df_q = pd.DataFrame(mdata1[:,:])
del mat1
del mdata1

###       Transforming input file of precipitation from .mat to DataFrame          ###
##     Sub-catchments/River segments
mat1 = loadmat("D:/working_directory/rainfall.mat")
mdata1 = mat1['...']
df_p = 0
df_p = pd.DataFrame(mdata1[:,:])
del mat1
del mdata1

###       Transforming input file of evapotranspiration from .mat to DataFrame          ###
##     Sub-catchments/River segments
mat1 = loadmat("D:/working_directory/evapo.mat")
mdata1 = mat1['...']
df_e = 0
df_e = pd.DataFrame(mdata1[:,:])
del mat1
del mdata1


###      Transforming input file of discharge from .xls to DataFrame        ###
##
#df_q = pd.read_excel(r'../discharge.xlsx', engine='openpyxl')

###      Transforming input file of precipitation from .xls to DataFrame        ###
##
#df_p = pd.read_excel(r'../rainfall.xlsx', engine='openpyxl')

###  Transforming input file of evapotranspiration from .xls to DataFrame ###
##
#df_e = pd.read_excel(r'../precipitation.xlsx', engine='openpyxl')


################################################################## 
###                SECTION 2 - CREATING DATE INFORMATION                    ###
##################################################################
           
### Creating DataFrama with date format e number of days of time step between entries
##
# Input intial date of the series of discharge

start = datetime.datetime.strptime("01/01/1980", "%d/%m/%Y").date() #TODO put start 
#data of flow file
numdays = len(df_q)
datelist = []
for x in range (0, numdays):
    datelist.append(start + datetime.timedelta(days = x))
    
# Creating a DataFrame with days, dates and discharges of a river segment
##
# Creating days list (number of days of the time step)
z = np.arange(len(df_q))
z = z+1

dias_contagem = list(datelist)


############################################################
###              SECTION 3 - CREATING EMPTY MATRIX FOR INDICES             ###
############################################################

######################################################
# Creates matrix to alocate calculated indices
z1 = np.arange(len(df_q.columns))
z1 = list(z1+1)


z2 = z1.copy()
del z2[1:(len(z1))]

z3 = z1.copy()
del z3[35:(len(z1))]

assinaturas = 0
assinaturas = pd.DataFrame(columns=['K90','IA','S','BFI',
                                    'STR_ela',"Q10_norm",'P_media','E_media'
                                    ],index=z1)

del z
del x
del numdays
del start

# Input year of the first entry in the series of discharges
yearlist = 0
yearlist = [datetime.datetime.strptime("1980", "%Y").year + index for index
            in range(len(z3))] ##TODO Input year of the first entry in the series of discharges


#######################################################     
###                    SECTION 4 - CALCULATING INDICES                  ###
#######################################################
    
########################### SEGMENT LOOP ###########################

before = datetime.datetime.now().isoformat(timespec='minutes')
before = pd.to_datetime(before) 
     
k=0
for k in range (0,(len(df_q.columns))):
    
        
    q_mini = list(df_q.iloc[:,k])
    p_mini = list(df_p.iloc[:,k])
    e_mini = list(df_e.iloc[:,k])
    
    matriz_q = 0
    matriz_q = pd.DataFrame({"dia": dias_contagem,"data": datelist,"vazao":q_mini})
    matriz_q['data'] = pd.to_datetime(matriz_q['data'])
    
    matriz_p = 0
    matriz_p = pd.DataFrame({"dia": dias_contagem,"data": datelist,"precip":p_mini})
    matriz_p['data'] = pd.to_datetime(matriz_p['data'])
    
    matriz_e = 0
    matriz_e = pd.DataFrame({"dia": dias_contagem,"data": datelist,"evapo":e_mini})
    matriz_e['data'] = pd.to_datetime(matriz_e['data'])
          
    print("Calculando trecho " + str(k+1) + " de " + str(len(df_q.columns)))
    
    
    ################## DISCHARGE RATE ##################
    # Calculating flow duration curve
    
    q_flow = matriz_q['vazao']
    q_flow = q_flow[np.logical_not(np.isnan(q_flow))]
    q_sort = pd.DataFrame(columns=['vazao'])
    q_sort['vazao'] = np.sort(q_flow)[::-1]
    q_sort['excedencia'] = np.arange(1.,len(q_sort)+1) / len(q_sort)
    q_sort['excedencia'] = q_sort['excedencia']*100

    #Calculatin discharge rate Q10
    q_10 = q_sort.iloc[(q_sort['excedencia']-10).abs().argsort()[:1]]
    
    #Calculatin discharge rate Q50
    q_50 = q_sort.iloc[(q_sort['excedencia']-50).abs().argsort()[:1]]

    #Calculatin discharge rate Q90
    q_90 = q_sort.iloc[(q_sort['excedencia']-90).abs().argsort()[:1]]

    
    ###########   1 - Coefficient of baseflow recession (k)  ###########
    ##

    # Calculatin ascension and descension rate
    
    print("   coeficiente de recessao da vazao de base (K)")     
    try:
        i = 0
        j = 0
                
        ano_g = 0
        ano_g = matriz_q.copy()
        ano_g['GD'] = 0
        ano_g['GA'] = 0
        ano_g['k'] = 0
        ano_g['k2'] = 0
        

        ano_g['GD'] = -1*(ano_g['vazao']-ano_g['vazao'].shift(-1))
        ano_g['GD'] = np.where(ano_g['GD']<0, ano_g['GD'], np.nan)
        ano_g['GA'] = -1*(ano_g['vazao']-ano_g['vazao'].shift(-1))
        ano_g['GA'] = np.where(ano_g['GA']>0, ano_g['GA'], np.nan)
        ano_g['k'] = np.where(ano_g['GD']<0,
                (-1/(np.log(ano_g['vazao'].shift(-1)/ano_g['vazao']))),
                              np.nan)
        
                            
        ano_g60 = ano_g.copy() 
        ano_g70 = ano_g.copy()
        ano_g80 = ano_g.copy()
        ano_g90 = ano_g.copy()


        #For discharges lower than Q90
        
        ano_g90['k'] = np.where(
            (ano_g90['vazao']>q_90.iloc[0,0]),
            np.nan, ano_g90['k'])
        
        ano_g90['k'] = np.where(
            (ano_g90['vazao']<0),
            np.nan, ano_g90['k'])
        
        array_g = pd.notnull(ano_g90['k'])

        maskg = [True, True, True, True, True]
        aux = np.zeros(len(array_g))
        
        for i in range(0, len(array_g)):
            chunk = array_g[i:i+5]
            if np.sum(chunk) == np.sum(maskg):
                aux[i:i+5] =+1
                
            else:
                pass
            
        aux = np.where((aux==0),np.nan, aux)  
        
        ano_g90['k2'] = ano_g90['k']*aux


    except IndexError:
        pass
    
    assinaturas.iloc[k,0] = abs(ano_g90['k2'].median())


    ###########            2 - Index of aridity (IA)                ###########
    ##
    ##
    # Calculating starting month of the hidrological year for minimums
    q = pd.DataFrame()
    q['data'] = matriz_q['data']
    q['data'] = q['data'].dt.year
    q= q.groupby('data').size()
    
    anos_min = 0
    anos_min = pd.DataFrame(index=list(yearlist),
                            columns = ['Data da vazao minima','Mes da vazao minima'])
    anos_min['Data da vazao minima'] = pd.to_datetime(anos_min['Data da vazao minima'])
      
    i = 0
    for i in range (0,len(yearlist)):
        ano1 = matriz_q.copy()
        ano1 = ano1.loc[(ano1['data'].dt.year==yearlist[i])]
        mini = ano1.copy()
        mini = mini.loc[(mini['vazao']==mini['vazao'].min())]
        
        if mini.size > 0:
            
            mini2 = mini.iloc[0] 
            anos_min.at[yearlist[i],'Data da vazao minima'] = mini2['data']
            anos_min['Mes da vazao minima'] = anos_min['Data da vazao minima'].dt.month
        else:
            pass
                    
        
    anos_min = anos_min.groupby('Mes da vazao minima').size().sort_values(ascending=False)
    mes_ano_min = anos_min.index[0]-5
    
    if mes_ano_min <= 0:
        mes_ano_min = mes_ano_min+12
    mes_anterior = mes_ano_min - 1
    if mes_anterior == 0:
        mes_anterior = 12

    # Deleting entries out of hydrological year range

    matriz_min = 0
    
    matriz_min = matriz_q.copy()
    
    i = 0
    for i in range (0,365):
        while matriz_min['data'][matriz_min.index[i]].month!=mes_ano_min:
            matriz_min = matriz_min.drop(index=[matriz_min.index[i]],axis=1)
        break
  
    i=0           
    matriz2_min = matriz_min.iloc[::-1]      
    for i in range (0,365):
        while matriz2_min['data'][matriz2_min.index[i]].month!=mes_anterior:
            matriz2_min = matriz2_min.drop(index=[matriz2_min.index[i]],axis=1)
        break
    matriz2_min = matriz2_min.iloc[::-1]

    yearlist2 = 0
    yearlist2 = matriz2_min['data'].dt.year
    yearlist2 = pd.DataFrame(yearlist2)
    yearlist2 = yearlist2.groupby('data').size()

    print("   indice de aridez (IA)")  

    e_media = 0
    e_media = pd.DataFrame(index = np.arange(len(yearlist2)), columns=['e_media',])
    p_media = 0
    p_media = pd.DataFrame(index = np.arange(len(yearlist2)), columns=['p_media',])
    q_media = 0
    q_media = pd.DataFrame(index = np.arange(len(yearlist2)), columns=['q_media',])
    

    try:
        i = 0
        for i in range (0,len(yearlist2)):
            mes = matriz2_min.loc[
                    matriz2_min['data'].dt.year == matriz2_min.iloc[0,1].year+i
                    ]
            mes = mes.loc[mes['data'].dt.month >= mes_ano_min]
            mes2 = matriz2_min.loc[
                    matriz2_min['data'].dt.year == matriz2_min.iloc[0,1].year+i+1
                    ]
            mes2 = mes2.loc[mes2['data'].dt.month < mes_ano_min]
            mes = mes.append(mes2)

            mask=(matriz_e['data']>=mes.iloc[0,1]) & (matriz_e['data']<=mes.iloc[len(mes)-1,1])

            evapo = matriz_e.loc[mask]

            e_media.iloc[i,0] = evapo['evapo'].mean()

            mask2=(matriz_p['data']>=mes.iloc[0,1]) & (matriz_p['data']<=mes.iloc[len(mes)-1,1])

            precip = matriz_p.loc[mask2]

            p_media.iloc[i,0] = precip['precip'].mean()
            
            mask3=(matriz_q['data']>=mes.iloc[0,1]) & (matriz_q['data']<=mes.iloc[len(mes)-1,1])

            vazao = matriz_q.loc[mask3]

            q_media.iloc[i,0] = vazao['vazao'].mean()
    
    except IndexError:
        pass


    # assinaturas.iloc[k,1] = (e_media['e_media'].mean())/(p_media['p_media'].mean())


    ###########  4    - Index of hydrogram assymetry (S)      ###########
    ##
    print("   indice de assimetria do hidrograma (S) ")

    S=1-(((1/ano_g['GD'].notnull().sum())*(ano_g['GD'].abs().sum()))/(
        (1/ano_g['GA'].notnull().sum())*(ano_g['GA'].abs().sum())))

    assinaturas.iloc[k,2] = S

    ###########            5    - Baseflow index (BFI)              ###########
    ##
    print("   baseflow index (BFI) ")
    
    BFI = q_90.iloc[0,0]/q_50.iloc[0,0]
    
    assinaturas.iloc[k,3] = BFI
    
    
    ###########              6    - Streamflow Elasticity                ###########
    ##
    
    print("   streamflow elasticity) ")
    
    st_el = (
        ((q_media-q_media.mean())['q_media']/
        (p_media-p_media.mean())['p_media']
        )*(
        p_media.mean()['p_media']/q_media.mean()['q_media'])
        ).median()
    
    assinaturas.iloc[k,4] = st_el    
    
    
    ###########              7    - Q10/Q50                         ###########
    ##
    
    print("   Q10/Q50 ")
    
    q10_norm = q_10.iloc[0,0]/q_50.iloc[0,0]
    
    assinaturas.iloc[k,5] = q10_norm
    
    
after = datetime.datetime.now().isoformat(timespec='minutes')
after = pd.to_datetime(after) 
tempo_execucao = after-before
print(tempo_execucao)

duration = 1000  # milliseconds
freq = 440  # Hz
winsound.Beep(freq, duration)


###########################################################
###            SECTION 5 - EXPORTING SIGNATURES RESULTS                     ###
###########################################################
   
################## Export DataFrame signature in .xls format ##################
writer = pd.ExcelWriter(r'D:/working_directory/results.xlsx')
assinaturas.to_excel(writer, sheet_name='Sheet1', index=False)
writer.save()

##########################        END         ###########################

