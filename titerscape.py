import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import sys, csv
import seaborn as sns
import os
import math
import re
import glob
import time
from itertools import repeat
from scipy import interpolate

def load_file(dict):      
    additional_list = ["POS1", "POS2"]
    pos = [3, 11]

    #generate column numbers
    col_id_num = []
    for i in range(1, 13):
        col_id_num.append("_C" + str(i))
    col_id_num = col_id_num * 2

    ID_list = [] 
    ID_index = []
    sub_data = {}
    for key, data in dict.items():
        #get list of id
        ID = data["ID"].dropna().tolist()
        ID = ID[len(ID) - ID[::-1].index("IgG") : ID.index("IgA")]
        acc = 0  #insert pos 1 & 2 to list
        for i in range(len(additional_list)):
            ID.insert(pos[i] + acc, additional_list[i])
            acc += 1
        ID = ["F" + str(key) + "_" + i for i in ID]     #Make Id unique
        ID = [i + j for i, j in zip(ID, col_id_num)]    #append column # to id
        ID_list.append(ID)         
        
        #get id index
        index = data.columns.get_loc("ID")
        data = data.iloc[:, index+1 : index+14].dropna()
        data.reset_index(drop=True, inplace=True)
        data.rename(columns={data.columns[0]: "Concentration"}, inplace = True) 
        sub_data[key] = data

    #flatten nested ID list to find duplicate IDs
    flatten = [item for sublist in ID_list for item in sublist]
    #replace special characters
    flatten = [re.sub("[:/\-()<>* ]","_",x) for x in flatten]

    duplicate = [idx for idx, item in enumerate(flatten) if item in flatten[:idx]]
    for i in range(len(duplicate)):
        flatten[duplicate[0]] = flatten[duplicate[0]] + "_" + str(i+1)
    sub_group = len(flatten) / len(input_data)
    ID_list = [flatten[i: i+int(sub_group)] for i in range(0, len(flatten), int(sub_group))]

    IgG, IgA, IgM = [], [], []
    for anti in range(len(data["Concentration"].unique())*2):  #mult by 2 bc of 2 sets
        IgG.append("IgG")
        IgA.append("IgA")
        IgM.append("IgM")
    anti_set = IgG + IgA + IgM
    anti_list = pd.DataFrame({"Antibody": anti_set})

    #create an antibody list for the graphs
    IgG_graph, IgA_graph, IgM_graph = [], [], []
    for anti in range(len(data["Concentration"].unique())):  #mult by 2 bc of 2 sets
        IgG_graph.append("IgG")
        IgA_graph.append("IgA")
        IgM_graph.append("IgM")
    anti_set_graph = IgG_graph + IgA_graph + IgM_graph
    anti_set_graph = pd.DataFrame({"Antibody": anti_set_graph})

    #twist the dataset; put plate 2 next to plate 1
    reindex = list(range(1,25)) 
    processed = {}
    for key, data in sub_data.items():
        data = pd.concat([anti_list, data ], axis=1)

        #create antibody subsets
        IgG_group = data[(data["Antibody"] == "IgG")]
        IgA_group = data[(data["Antibody"] == "IgA")]
        IgM_group = data[(data["Antibody"] == "IgM")]

        IgG_unique = IgG_group.drop_duplicates("Concentration")
        IgG_unique.reset_index(drop=True, inplace=True)
        IgG_duplicate = IgG_group[IgG_group.duplicated("Concentration")]
        IgG_duplicate = IgG_duplicate.drop(IgG_duplicate.columns[0:2], axis=1)
        IgG_duplicate.reset_index(drop=True, inplace=True)
        IgG_list = pd.concat([IgG_unique, IgG_duplicate], axis=1, sort=False)

        IgA_unique = IgA_group.drop_duplicates("Concentration")
        IgA_unique.reset_index(drop=True, inplace=True)
        IgA_duplicate = IgA_group[IgA_group.duplicated("Concentration")]  
        IgA_duplicate = IgA_duplicate.drop(IgA_duplicate.columns[0:2], axis=1)
        IgA_duplicate.reset_index(drop=True, inplace=True)
        IgA_list = pd.concat([IgA_unique, IgA_duplicate], axis=1, sort=False)

        IgM_unique = IgM_group.drop_duplicates("Concentration")
        IgM_unique.reset_index(drop=True, inplace=True)
        IgM_duplicate = IgM_group[IgM_group.duplicated("Concentration")]
        IgM_duplicate = IgM_duplicate.drop(IgM_duplicate.columns[0:2], axis=1)
        IgM_duplicate.reset_index(drop=True, inplace=True)
        IgM_list = pd.concat([IgM_unique, IgM_duplicate], axis=1, sort=False)

        result = pd.concat([IgG_list, IgA_list, IgM_list], sort=False)

        result.columns = result.columns[:2].tolist() + reindex
        result.reset_index(drop=True, inplace=True)
        result = result[result.columns[2:]]
        result = result.T
        result.reset_index(drop=True, inplace=True)

        processed[key] = result

    #generate concentration list
    concentration = []
    for i in range(7):
        if i == 0:
            c = (1/50)
            concentration.append(c)
        else:
            c = c/3
            concentration.append(c)
    concentration.append(0)

    #generate double headers
    header = pd.MultiIndex.from_product([["IgG", "IgA", "IgM"],concentration])

    #generate column numbers
    col_num = []
    for i in range(1, 13):
        col_num.append(i)
    col_num = col_num * 2
    col_num = pd.DataFrame({"Column": col_num})

    output = {}
    input_name = []
    for key, data in processed.items():
        data.columns = header    
        ID = pd.DataFrame({"Sample ID": ID_list[key-1]})

        input_name.append(re.split('V-2|\+', str(input_data[key-1]))[1]) #retrive file names

        result = pd.concat([col_num,  ID, data], axis=1)
        output[key] = result

    #graph_table = output.copy()

    #append all dataframes into one
    result_data = pd.DataFrame()

    for key, sub_df in output.items():
        result_data = result_data.append(sub_df, ignore_index = False)

    result_data.reset_index(drop=True, inplace=True)

    n = len(result_data) / len(input_name)
    name = [item for item in input_name for i in range(int(n))]
    file = pd.DataFrame({"File": name}) 
    file.reset_index(drop=True, inplace=True)
    result_data = pd.concat([file, result_data], axis = 1)

    #create a duplicate table
    graph_table = result_data.copy()

    concent_set = concentration * 3
    concent_list = [round(num, 4) for num in concent_set]
    concent_list = [str(i) for i in concent_list]
    concent_list.insert(0, "Sample ID")

    #making graphs
    graph_table = graph_table.T
    new_header = graph_table.iloc[2]
    graph_table = graph_table[3:]
    graph_table.columns = new_header

    concent_set = pd.DataFrame({"Log(DF)": concent_set})

    #reset index prior to concatenation 
    concent_set.reset_index(drop=True, inplace=True)
    anti_set_graph.reset_index(drop=True, inplace=True)
    graph_table.reset_index(drop=True, inplace=True)
    graph_table = pd.concat([concent_set, anti_set_graph, graph_table], axis = 1).reset_index(drop=True)
    graph_table_sig = graph_table.copy() #duplicate dataset

    if graph_table["Log(DF)"].any() != 0:
        graph_table["Log(DF)"] = np.log10(graph_table["Log(DF)"])
        graph_table.replace([np.inf, -np.inf], np.nan, inplace=True)  #cannot take log of 0; replace it with NA
        graph_table = graph_table[graph_table["Log(DF)"].notna()]  #drop rows with NA   

    if graph_table_sig["Log(DF)"].any() != 0:
        graph_table_sig["Log(DF)"] = np.log10(graph_table_sig["Log(DF)"])
        graph_table_sig.replace([np.inf, -np.inf], 0, inplace=True)  #cannot take log of 0; replace it with NA      

    #convert column values to numeric 
    graph_table.iloc[:, 2:] = graph_table.iloc[:, 2:].astype(float)
    #create antibody subsets
    IgG_group = graph_table[(graph_table["Antibody"] == "IgG")]
    IgA_group = graph_table[(graph_table["Antibody"] == "IgA")]
    IgM_group = graph_table[(graph_table["Antibody"] == "IgM")]

    return result_data, graph_table, IgG_group, IgA_group, IgM_group, graph_table_sig

def output(result_data, graph_table_sig, IgG_cutoff, IgA_cutoff, IgM_cutoff):
    #convert column values to numeric 
    graph_table_sig = graph_table_sig.loc[graph_table_sig['Log(DF)'] != 0]
    graph_table_sig.iloc[:, 2:] = graph_table_sig.iloc[:, 2:].astype(float)

    IgG_group_sig = graph_table_sig[(graph_table_sig["Antibody"] == "IgG")]
    IgA_group_sig = graph_table_sig[(graph_table_sig["Antibody"] == "IgA")]
    IgM_group_sig = graph_table_sig[(graph_table_sig["Antibody"] == "IgM")]

    IgG_endpoint, IgA_endpoint, IgM_endpoint = [], [], []
    for col in IgG_group_sig.columns[2:len(IgG_group_sig.columns)]:  
        x = list(IgG_group_sig["Log(DF)"])
        y = list(IgG_group_sig[col])
        
        if IgG_cutoff > max(y):
        	IgG_endpoint.append(0)
        else: 
        	ac = interpolate.interp1d(y, x, fill_value='extrapolate')
        	IgG_endpoint.append(ac(IgG_cutoff) - (ac(IgG_cutoff)* 0.01259640486))
        
    IgG_act_val = IgG_endpoint[:]
    for i in range(len(IgG_act_val)):
        IgG_act_val[i] = 1 / (10**IgG_act_val[i])               #math.pow(10, IgG_act_val[i])
            
    for col in IgA_group_sig.columns[2:len(IgA_group_sig.columns)]:  
        x = list(IgA_group_sig["Log(DF)"])
        y = list(IgA_group_sig[col])

        if IgA_cutoff > max(y):
        	IgA_endpoint.append(0)
        else:
        	ac = interpolate.interp1d(y, x, fill_value='extrapolate')
        	IgA_endpoint.append(ac(IgA_cutoff) - (ac(IgA_cutoff)* 0.01259640486))  #- (ac(IgA_cutoff)* 0.01259640486)
        
    IgA_act_val = IgA_endpoint[:]
    for i in range(len(IgA_act_val)):
        IgA_act_val[i] = 1 / (10**IgA_act_val[i])                 #math.pow(10, IgA_act_val[i])
        
    for col in IgM_group_sig.columns[2:len(IgM_group_sig.columns)]:  
        x = list(IgM_group_sig["Log(DF)"])
        y = list(IgM_group_sig[col])
        
        if IgM_cutoff > max(y):
        	IgM_endpoint.append(0)
        else:
        	ac = interpolate.interp1d(y, x, fill_value='extrapolate')
        	IgM_endpoint.append(ac(IgM_cutoff) - (ac(IgM_cutoff)* 0.01259640486)) 
        
    IgM_act_val = IgM_endpoint[:]
    for i in range(len(IgM_act_val)):
        IgM_act_val[i] = 1 / (10**IgM_act_val[i])                        #math.pow(10, IgM_act_val[i])

    IgG_endpoint = pd.DataFrame({"IgG Log Endpoint": IgG_endpoint})
    IgG_act_val = pd.DataFrame({"IgG Endpoint": IgG_act_val})
    IgA_endpoint = pd.DataFrame({"IgA Log Endpoint": IgA_endpoint})
    IgA_act_val = pd.DataFrame({"IgA Endpoint": IgA_act_val})
    IgM_endpoint = pd.DataFrame({"IgM Log Endpoint": IgM_endpoint})
    IgM_act_val = pd.DataFrame({"IgM Endpoint": IgM_act_val})
    result_data = pd.concat([result_data, IgG_endpoint, IgG_act_val, IgA_endpoint, IgA_act_val, IgM_endpoint, IgM_act_val], axis=1)

    return result_data

# x0 = the minimum value that can be obtained (i.e. what happens at 0 dose)
# y0 = the maximum value that can be obtained (i.e. what happens at infinite dose)
# c = the point of inflection (i.e. the point on the S shaped curve halfway between x0 and y0)
# k = Hill’s slope of the curve (i.e. this is related to the steepness of the curve at point c).

#sigmoidal functions
def sigmoid(p,x):
    x0,y0,c,k=p
    y = c / (1 + np.exp(-k*(x-x0))) + y0
    return y

def residuals(p,x,y):
    return y - sigmoid(p,x)

#use resize if normalized is needed
def resize(arr,lower=0.0,upper=1.0):    #lower=0.0, upper=1.0
    arr=arr.copy()
    if lower>upper: lower,upper=upper,lower
    arr -= arr.min()
    arr *= (upper-lower)/arr.max()
    arr += lower
    return arr

#calculate endpoint titer
class AxisCam:
    def __init__(self, x=None, y=None):
        self.x = x if x else []
        self.y = y if y else []

        if len(self.x):
            self.xMin = min(self.x)
            self.xMax = max(self.x)
        else:
            self.xMin = None
            self.xMax = None

        if len(self.y):
            self.yMin = min(self.y)
            self.yMax = max(self.y)
        else:
            self.yMin = None
            self.yMax = None

        self._interpolX, self._interpolY = self.setInterpolator()

    def setInterpolator(self, interpolator=interpolate.interp1d):
        """
        Define the interpolator to use to approximate the axis cam positions
        :param interpolator: interpolator function to use, default is scipy.interpolate.interp1d
        :return: a tuple with the interpolator functions for x and y values
        """
        if len(self.x) <= 0 or len(self.y) <= 0:
            return None, None
        with np.errstate(divide='ignore', invalid='ignore'):  # silent the warnings caused by the interpolator
            self._interpolX = interpolator(self.y, self.x)  # x = f(y)
            self._interpolY = interpolator(self.x, self.y)  # y = f(x)
        return self._interpolX, self._interpolY

    def getX(self, yValue):
        """
        Return x-value corresponding to a y-value using the interpolator
        :param yValue: y-value we want to know the corresponding x-value
        :return: x-value corresponding to the given y-value
        """
        if yValue < self.yMin:
            #raise ValueError("value should be greater than the minimum y-value")
            return 0
        elif yValue > self.yMax:
            #raise ValueError("value should be lesser than the maximum y-value")
            return 0
        return float(self._interpolX(yValue))

    def getY(self, value):
        """
        Return a y-value corresponding to a x-value using the interpolator
        :param value: x-value we want to know the corresponding y-value
        :return: the y-value corresponding to the given x-value
        """
        if value < self.xMin:
            raise ValueError("value should be greater than the minimum x-value")
        elif value > self.xMax:
            raise ValueError("value should be lesser than the maximum x-value")
        return float(self._interpolY(value))

#Main Method
if __name__ == '__main__':

    if len(sys.argv) < 5:
        print("Expected input format: python covid_v3.py <output file name><IgG cutoff><IgA cutoff><IgM cutoff>")
        sys.exit("Please try again.")
    else:
        output_file = sys.argv[1]
        IgG_cutoff = float(sys.argv[2])
        IgA_cutoff = float(sys.argv[3])
        IgM_cutoff = float(sys.argv[4])

    start = time.time()

    input_data = glob.glob("*.xlsx")
    dict = {}
    for i in range(len(input_data)):
        dict[i+1] = pd.read_excel(input_data[i]) 

    result_data, graph_table, IgG_group, IgA_group, IgM_group, graph_table_sig = load_file(dict) 
    result_data = output(result_data, graph_table_sig, IgG_cutoff, IgA_cutoff, IgM_cutoff)
    result_data.to_csv(output_file, index=False)

    #create directory for figs
    script_dir = os.path.dirname(__file__)
    fig_dir = os.path.join(script_dir, 'fig/')
    #for sigmoidal curves
    fig_IgG_dir = os.path.join(script_dir, 'fig/IgG/') 
    fig_IgA_dir = os.path.join(script_dir, 'fig/IgA/')
    fig_IgM_dir = os.path.join(script_dir, 'fig/IgM/')

    #if directory path does not exist, make one
    if not os.path.isdir(fig_dir):
        os.makedirs(fig_dir)
    if not os.path.isdir(fig_IgG_dir):      
        os.makedirs(fig_IgG_dir)
    if not os.path.isdir(fig_IgA_dir):      
        os.makedirs(fig_IgA_dir)
    if not os.path.isdir(fig_IgM_dir):      
        os.makedirs(fig_IgM_dir)

    #make main figures
    for col in graph_table.columns[2:len(graph_table.columns)]:   
        sns.lineplot(x="Log(DF)", y=col, hue="Antibody", data=graph_table, marker="o").set_title(col)  
       # sns.despine(offset=10, trim=True)
        plt.legend(loc='upper right')
        plt.ylabel('OD490')
        plt.xlim(-5, -1.5)
        plt.gca().invert_xaxis() #reverse x-axis
        plt.savefig(fig_dir + col + ".png")
        plt.clf()
    plt.close()

    #make sigmoidal curves for each antibody
    #IgG
    for col in IgG_group.columns[2:len(IgG_group.columns)]:  
        x = IgG_group["Log(DF)"]
        y = IgG_group[col]
        
        p_guess=(np.median(x),np.median(y),1.0,1.0)
        p, cov, infodict, mesg, ier = scipy.optimize.leastsq(
            residuals,p_guess,args=(x,y),full_output=1) 

        x0,y0,c,k=p
        #print('''\
        #x0 = {x0}
        #y0 = {y0}
        #c = {c}
        #k = {k}
        #'''.format(x0=x0,y0=y0,c=c,k=k))
        xp = np.linspace(-4.6, -1.5, 1500)
        pxp=sigmoid(p,xp)

        # Plot the results
        plt.plot(x, y, '.', xp, pxp, '-')
        plt.gca().invert_xaxis() #reverse x-axis
        #plt.axhline(y=0.5, color='r', linestyle='-')
        plt.suptitle(col)
        plt.xlabel('Log(DF)')
        plt.ylabel('OD490')  #,rotation='horizontal'
        plt.savefig(fig_IgG_dir + col + ".png")
        plt.clf()
    plt.close()
    #IgA
    for col in IgA_group.columns[2:len(IgA_group.columns)]:  
        x = IgA_group["Log(DF)"]
        y = IgA_group[col]
        
        p_guess=(np.median(x),np.median(y),1.0,1.0)
        p, cov, infodict, mesg, ier = scipy.optimize.leastsq(
            residuals,p_guess,args=(x,y),full_output=1) 

        x0,y0,c,k=p
        xp = np.linspace(-4.6, -1.5, 1500)
        pxp=sigmoid(p,xp)

        # Plot the results
        plt.plot(x, y, '.', xp, pxp, '-')
        plt.gca().invert_xaxis() #reverse x-axis
        #plt.axhline(y=0.5, color='r', linestyle='-')
        plt.suptitle(col)
        plt.xlabel('Log(DF)')
        plt.ylabel('OD490')  #,rotation='horizontal'
        plt.savefig(fig_IgA_dir + col + ".png")
        plt.clf()
    plt.close()
    #IgM
    for col in IgM_group.columns[2:len(IgM_group.columns)]:  
        x = IgM_group["Log(DF)"]
        y = IgM_group[col]
        
        p_guess=(np.median(x),np.median(y),1.0,1.0)
        p, cov, infodict, mesg, ier = scipy.optimize.leastsq(
            residuals,p_guess,args=(x,y),full_output=1) 

        x0,y0,c,k=p
        xp = np.linspace(-4.6, -1.5, 1500)
        pxp=sigmoid(p,xp)

        # Plot the results
        plt.plot(x, y, '.', xp, pxp, '-')
        plt.gca().invert_xaxis() #reverse x-axis
        #plt.axhline(y=0.5, color='r', linestyle='-')
        plt.suptitle(col)
        plt.xlabel('Log(DF)')
        plt.ylabel('OD490')  #,rotation='horizontal'
        plt.savefig(fig_IgM_dir + col + ".png")
        plt.clf()
    plt.close()

    end = time.time()
    sample_count = 24
    print("The runtime for " + str(len(input_data)) + " file(s) (" + str(sample_count * len(input_data)) +  " samples) is: " + str(end-start) + " sec")