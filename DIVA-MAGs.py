#!/usr/bin/env python

#DIVA-MAGs

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
import argparse
import os
from sys import exit

#argparse inputs, only --files is required
parser = argparse.ArgumentParser(description='DIVA-MAGs: Analysis of distribution of variants across MAGs')
parser.add_argument('--files', type=argparse.FileType('r'), help=".txt file with paths of the input files: scaffold_info.tsv, scaffold_to_bin.stb, SNVs.tsv files", required = True)
parser.add_argument("--coverage_info", "--c", action="store_true", help='Ask for information to decide threshold for coverage, by deafult the program use the first quartile of the coverage distribution')
parser.add_argument("--window_info", "--w", action="store_true", help='Ask for MAGs information to decide the window size and step, if not given then the program manually create window size and step for each MAG')
parser.add_argument("--printed_position", "--p", action = "store", type = int, default = 20, help='Range of position to print in the final output, default 20')
parser.add_argument("--min_variants_MAG", "--min-MAG", action = "store", type = int, default = 20, help="Minimum number of variants needed to perform MAG analysis, otherwise discarded, by default 20")
parser.add_argument("--min_variants_scaffold", "--min-scaffold", action = "store", type = int, default = 1, help="Minimum number of variants needed to perform scaffold analysis, otherwise discarded. Greater than 1, by default 1")
args = parser.parse_args()

#read the input files 
def deal_with_inputs(path_scaffold_info,path_scaffold_to_bin,path_SNVs):
    scaffold_info = pd.read_csv(path_scaffold_info, sep = '\t', index_col = "scaffold")
    scaffold_to_bin = pd.read_csv(path_scaffold_to_bin, sep = '\t', header = None, names = ["scaffold", "MAG"], index_col = "scaffold")
    #scaffold2bin creation: scaffold, length, MAG as columns
    scaffold2bin = scaffold_to_bin.join(scaffold_info["length"])
    SNVs = pd.read_csv(path_SNVs, sep = "\t", index_col = "scaffold")
    return scaffold2bin,SNVs
    
paths_inputs = []
for path in args.files:
    paths_inputs.append(path.strip())
sample_name,path_scaffold_info,path_scaffold_to_bin,path_SNVs = paths_inputs

files= deal_with_inputs(path_scaffold_info,path_scaffold_to_bin,path_SNVs)
scaffold2bin, SNVs = files


#asking for directory for storing the results
path_res = input("Give me a path to save the output: ")
print("\t")
#check existence, create new folder for the results
if os.path.exists(path_res) == False:
    print('Error: not existing path')
    exit('Game over!')
else:
    path = os.path.join(path_res, sample_name)
    os.makedirs(path, exist_ok = True)  
    os.chdir(path)

#giving information on coverage distribution if asked by user 
if args.coverage_info:
    print('Printing information of the coverage distribution to help to decide filtering threshold \n')
    print(pd.DataFrame(SNVs["position_coverage"].describe()[1:]))
    print("\t")
    #asking for coverage threshold (int) to filter scaffolds
    threshold = input('Please give us the position coverage threshold for filtering ')
    try:
        int(threshold)
    except ValueError:
        print('Error: threshold for position coverage must be an integer')
        exit('Game over!')
else:
    threshold = int(SNVs["position_coverage"].describe()[4])

#filtering scaffolds for coverage  
if threshold != 0:
    print("Filtering for the coverage, threshold:", threshold, '\n')
    filter_SNVs = SNVs[SNVs["position_coverage"] > int(threshold)]
    filtered_scaffold2bin = scaffold2bin
    filtered = {'MAG':[], 'scaffold':[]}
    for scaffold in scaffold2bin.index.values:
        if scaffold not in filter_SNVs.index.values:
            filtered_scaffold2bin = filtered_scaffold2bin.drop([scaffold])
            filtered['MAG'].append(scaffold2bin.loc[scaffold, 'MAG'])
            filtered['scaffold'].append(scaffold)
    scaffold2bin = filtered_scaffold2bin
    print("Number of scaffold after filtering is ", len(scaffold2bin))
    print('List of filtered scaffolds for coverage available on filtered_scaffolds.csv file')
    df_filtered = pd.DataFrame(data=filtered)
    df_filtered.to_csv("filtered_scaffolds.csv")
    filter_SNVs_MAG = filter_SNVs.join(scaffold2bin["MAG"])
else:
    filter_SNVs_MAG = SNVs.join(scaffold2bin["MAG"])

#creating sub-dataframes of SNVs and scaffold2bin for MAGs
#dataframes stored in two dictionary, keys are MAG names
dic_filter_SNVs_MAG = {}
for MAG in filter_SNVs_MAG["MAG"].unique():
    dic_filter_SNVs_MAG[MAG] = filter_SNVs_MAG[filter_SNVs_MAG["MAG"] == MAG]   
    
dic_scaffold2bin = {}
for MAG in scaffold2bin["MAG"].unique():
    dic_scaffold2bin[MAG] = scaffold2bin[scaffold2bin["MAG"] == MAG].sort_values(by = ["length"])   
if threshold != 0:
    print("Number of MAGs after filtering is ", len(dic_scaffold2bin.keys()), '\n')

#giving information on MAGs structure to choose window sizes and step if asked by the user
#creating dataframe with MAG, scaffold length informations, variants informations
#asking for file with windows size and steps decided by user
def MAG_info():
    MAG_names = list(dic_filter_SNVs_MAG.keys())
    scaffold_number = []    
    Avg_scaffold_len = []
    Min_scaffold_len = []
    Max_scaffold_len = []
    for MAG in dic_filter_SNVs_MAG.keys():
        scaffold_number.append(len(dic_scaffold2bin[MAG]))
        avg = np.mean(list(dic_scaffold2bin[MAG].loc[:,'length'].values))
        minimum = min(list(dic_scaffold2bin[MAG].loc[:,'length'].values))
        maximum = max(list(dic_scaffold2bin[MAG].loc[:,'length'].values))
        Avg_scaffold_len.append(avg)
        Min_scaffold_len.append(minimum)
        Max_scaffold_len.append(maximum)

    Avg_variants = []
    Min_variants = []
    Max_variants = []
    for MAG in dic_filter_SNVs_MAG.keys():
        length_data = []
        for scaffold in list(dic_filter_SNVs_MAG[MAG].index.unique()):
            length = len(dic_filter_SNVs_MAG[MAG].loc[scaffold,:])
            length_data.append(length)
        Avg_variants.append(np.mean(length_data))
        Min_variants.append(min(length_data))
        Max_variants.append(max(length_data))

    d = {'MAG': MAG_names, 'Scaffold number': scaffold_number, 'Avg scaffold length': Avg_scaffold_len, 'Min scaffold length': Min_scaffold_len, 'Max scaffold length': Max_scaffold_len, 
    'Avg variants number':Avg_variants, 'Min variants number': Min_variants, 'Max variants number': Max_variants}
    MAG_info = pd.DataFrame(data=d)
    MAG_info = MAG_info.set_index('MAG')
    MAG_info.to_csv("MAG_info.csv")
      
if args.window_info:
    MAG_info()
    #file must contain a row for each MAG
    #each row must contain window and step values tab-separated
    filepath = input('Enter your file with window and step here: ')
    print('\t')
    dic_windows = {}
    dic_steps = {} 
    count = 0
    #reading information from the file if provided correctly
    if os.path.exists(filepath) == False:
        print("Error: file don't exist")
        exit('Game over!')
    with open(filepath) as fp:    
        for MAG,line in zip(dic_filter_SNVs_MAG.keys(),fp):
            count += 1
            try:
                int((line.split(' ')[0]))
                int((line.split(' ')[1]))
            except ValueError:
                print('Error: window and step must two be integers')
                exit('Game over!')
            dic_windows[MAG] = line.split(' ')[0]
            dic_steps[MAG] = line.split(' ')[1]   
    if count != len(dic_filter_SNVs_MAG.keys()):
        print('Error: wrong number of inputs provided, window size and step per MAG are needed')
        exit('Game over!') 

#computing windows size per MAG if user doesn't decide it manually
def window_size(genome):
    window = 0
    n = 0
    for scaffold in list(genome.index.unique()):
        length = int(scaffold2bin.loc[scaffold, "length"])
        variant_count = int(genome.loc[scaffold, :].shape[0])
        window += ((length // variant_count) * 50)
        n += 1
    window = window // n
    return window

#computing number of variants per windows in scaffold 
def variant_window(scaffold):
    start_pos = 0
    last_pos = int(scaffold2bin.loc[scaffold, "length"])
    #dictionary with initial window position as key and variants number as value
    dic_windows = {}
    while (last_pos - start_pos) >= window:
        dic_windows[start_pos] = len(scaffold_dataframe[(scaffold_dataframe["position"] >= start_pos) & (scaffold_dataframe["position"] <= (start_pos + window))])
        start_pos += skip
    return dic_windows

#computing average of windows variants per position
def variant_position(scaffold):
    dic_variant_position = {}
    dic_windows = variant_window(scaffold)
    for i in range(int(scaffold2bin.loc[scaffold, "length"])):
        #check number of windows that span a position and store avg of values in dictionary
        count = 0
        number_window = 0
        #print only positions required by user, default 20
        if i % args.printed_position == 0:
            for key,value in dic_windows.items():
                if i >= key and i <= (key + window):
                    number_window += 1
                    count += value
        if number_window != 0:
            avg = count / number_window
            dic_variant_position[i] = [round(avg,2)]
    #return dataframe: position, average and scaffold
    final = pd.DataFrame(dic_variant_position)
    final = final.transpose()
    final = final.rename(columns = {0 : "average"})
    final["scaffold"] = scaffold
    final = final.reset_index()
    final = final.set_index("scaffold")
    final = final.rename(columns = {"index" : "position"})
    return final

#computing average of windows variants per position for each scaffold in each MAG
def variant_distribution():
    discarded = {'MAG':[],'scaffold':[]}
    windows = {'MAG':[],'window':[]}
    #check if MAG has enough variants
    for MAG in list(dic_filter_SNVs_MAG.keys()):
        print("Analysing",MAG, "...")
        if dic_filter_SNVs_MAG[MAG].shape[0] <= args.min_variants_MAG:
            print("Discarding ", MAG, ": insufficient number of variants \n")
            continue
        #assign windows and step value according to user input
        global genome
        genome = dic_filter_SNVs_MAG[MAG]
        global window
        if args.window_info:
            window = int(dic_windows[MAG])
        else:
            window = window_size(genome)
            print('Calculating window size and step values')
        windows['MAG'].append(MAG)
        windows['window'].append(window)
        global skip
        if args.window_info:
            skip = int(dic_steps[MAG])
        else:
            #if not provided by user, step is window // 10
            skip = window // 10
        print('Using window', window, 'and step', skip, 'on', MAG)
        #creating dataframe to store values 
        a = {"scaffold" : [], "position" : [], "average" : []}
        df = pd.DataFrame(a)
        df = df.set_index("scaffold")
        print('Scaffold analysis on', MAG)
        print('Discarding scaffold whose length is shorter than window size...')
        print('Discarding scaffold with low number of variants...')
        #checking if each scaffold can be analysed and call variant_position(scaffold)
        for scaffold in list(genome.index.unique()):
            if int(scaffold2bin.loc[scaffold, "length"]) >= window:
                global scaffold_dataframe
                scaffold_dataframe = genome.loc[scaffold, :]
                if type(scaffold_dataframe) != pd.DataFrame:
                    discarded['MAG'].append(MAG)
                    discarded['scaffold'].append(scaffold)
                    continue
                if args.min_variants_scaffold:
                    if len(scaffold_dataframe) < args.min_variants_scaffold:
                        discarded['MAG'].append(MAG)
                        discarded['scaffold'].append(scaffold)
                        continue
                df = df.append(variant_position(scaffold))
            else:
                discarded['MAG'].append(MAG)
                discarded['scaffold'].append(scaffold)
        #create .csv file for each MAG
        name = MAG
        print('Creating .csv file with variant informations')
        df.to_csv(f"{name}.csv")
        l = 0
        unique_scaffolds = df.index.unique()
        #create plot for each MAG
        print("Computing the MAG plot \n")
        for us in unique_scaffolds:
            S = df[df.index == us]
            plt.rcParams["figure.figsize"] = (30, 10)
            plt.plot(S['position'].values + l ,S['average'],label=str(us))
            plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
            plt.axvline(x = S['position'].values[0] + l)
            l += scaffold2bin.loc[us, "length"]
        plt.xlabel("position")
        plt.ylabel("variation average")
        plt.title(MAG)
        plt.savefig(f"{MAG}_plot.pdf")
        plt.close()
    #store window and discarded scaffold information in .csv files    
    discarded_df = pd.DataFrame(data=discarded)
    discarded_df.to_csv("Discarded_scaffolds.csv")
    windows_df = pd.DataFrame(data=windows)
    windows_df.to_csv("Windows.csv")
    print('Used windows size available in Windows.csv file')
    print('Information on discarded scaffold available in discarded.csv file')
    print(len(discarded_df), 'discarded scaffold \n') 
    return    

#call the function and store outputs in the created folder 
variant_distribution()

print('Thanks for being patient, goodbye!')


