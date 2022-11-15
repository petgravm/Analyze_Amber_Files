import gzip
import os
import time
import datetime
import pathlib
import subprocess
from stat import S_ISREG, ST_CTIME, ST_MODE
import sys
import re
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pkg_resources

#path= "C:/Users/mayap/Downloads/6_Runs/Hassan/"

           #CHOICE 1
def ChoiceofOutput(choice):
    path=input("Directory of output files:")
    Type=input("Molecular Dynamic Program being used, ex) AMBER, NAMD, etc:")
    if Type=='AMBER':
        Recent_File_AMBER(path, Type, choice)
        
    #elif: ......
 
    #Collects all files in the Directory Choice 2
def Recent_File_AMBER(path,Type, choice):
    x=""
    for filename in os.listdir(path): # os.listdir(args[path]) # and args[type]==AMBER
        if filename.endswith('.out') and filename !=x and Type=='AMBER' and int(choice)==2: #only go to file if it ends in '.out'--specifices for output files 
            entries = (os.path.join(path,filename)) #Get the directory
            fname=pathlib.Path(entries)
            print(fname)
    return(fname)


    #Choice 2: Individually inputting in the files
def GetFiles(choice):
    path=input("Directory of output files:")
    Type=input("Molecular Dynamic Program being used, ex) AMBER, NAMD, etc:")
    files=input("List the file names with extension:")
    if Type=='AMBER':
        Give_File_Amber(path, Type, choice,files)
    #elif ....
    
    
def Give_File_Amber(path, Type, choice, files):
    list_of_files= files.split(' ')
    num_files=len(list_of_files) #remaining files left to read
    num_of_files=len(list_of_files) #the total number of files that are read in
    #num_files=0
    #print(files.split(' '))
    #1) go into the directory and find each file in list_of_files
    
    data_list = []
    
    for file in list_of_files:
        if file in os.listdir(path):
            num_files=num_files-1
            entries = (os.path.join(path,file))
            print(entries)
            #num_files+=1
           # print(num_files)
            
            data_list.append(parse_amber(entries, num_files, num_of_files))
    
    final_dict = combine_dicts(data_list)
    
    Generate_Graph(final_dict)
    #print(final_dict['time'])
    
    #2) read through each file and call parse_amber()
    
def combine_dicts(list_of_dicts):
    
    for k, v in list_of_dicts[0].items():
        for i in range(1, len(list_of_dicts)):
            list_of_dicts[0][k].extend(list_of_dicts[i][k])
            
    return list_of_dicts[0]
    
    
def parse_amber(entries,num_files,num_of_files):
    '''Parse AMBER output file (.out), and return data as dict.'''
    #filepath='C:/Users/mayap/Downloads/6_Runs/Hassan/prod30.out'
    filepath=entries
    num=num_files
    total_files=num_of_files
    combined={}
    #num=num-1
    #print(num)
    data = dict(
        nstep=[], time=[], temp=[], 
        press=[], etot=[], ektot=[], 
        eptot=[], bond=[], angle=[], 
        dihed=[], ub=[], imp=[], cmap=[], 
        nb=[], eel=[], vdwaals=[], 
        eelec=[], ehbond=[], restraint=[], 
        ekcmt=[], virial=[], volume=[], 
        density=[], 
    )
    try:
        with open(filepath, 'r') as f:
            header = True
            for line in f.readlines():
                split_line = re.split('\s+', line)
                if header and 'NSTEP' in split_line:
                    header = False
                if not header:
                    if 'NSTEP' in split_line:
                        data['nstep'].append(split_line[3])
                        data['time'].append(split_line[6])
                        data['temp'].append(split_line[9])
                        data['press'].append(split_line[12])
                    elif 'Etot' in split_line:
                        data['etot'].append(split_line[3])
                        data['ektot'].append(split_line[6])
                        data['eptot'].append(split_line[9])
                    elif 'BOND' in split_line:
                        data['bond'].append(split_line[3])
                        data['angle'].append(split_line[6])
                        data['dihed'].append(split_line[9])
                    elif 'UB' in split_line:
                        data['ub'].append(split_line[3])
                        data['imp'].append(split_line[6])
                        data['cmap'].append(split_line[9])
                    elif '1-4' in split_line:
                        data['nb'].append(split_line[4])
                        data['eel'].append(split_line[8])
                        data['vdwaals'].append(split_line[11])
                    elif 'EELEC' in split_line:
                        data['eelec'].append(split_line[3])
                        data['ehbond'].append(split_line[6])
                        data['restraint'].append(split_line[9])
                    elif 'EKCMT' in split_line:
                        data['ekcmt'].append(split_line[3])
                        data['virial'].append(split_line[6])
                        data['volume'].append(split_line[9])
                    elif 'Density' in split_line:
                        data['density'].append(split_line[3])
                    elif 'A' and 'V' and 'E'and 'R' and 'A' and 'G' and 'E' and 'S' in split_line:
                        break
                        
    except FileNotFoundError:
        if filepath.contains('/'):
            filename = filepath.split('/')[-1]
        else:
            filename = 'provided input file'
        print(f'Warning: Could not find {filename}!')
        
    #Checks to see if parameter is present in output file
    #If it is not present, it removes the parameter and 
    #returns the dictonary without that parameter 
    
    #print(data['nstep'])
    
    for k,v in list(data.items()): #looks through the key('nstep') and values('600.2') in the data dictonary that holds each parameter
        if v==[]:  #v=[] checks if the value is empty or if the parameter doesn't have any values, and thus is not present within the file
            data.pop(k) #remove the parameter from the dictonary if it does not exist
            print('removing parameter: '+str(k))    
        else:
            continue       #if the parameter is not empty, continue...
            
    updated_dict=dict(data)  #gets the updated dictonary (with the removed empty parameters)
    print(len(data))
    return updated_dict
    #MergeDict(data,total_files,num)
    #if num_files==0:
        #Generate_Graph_1(updated_dict) 


def Generate_Graph(data):
    #make the graph
    files=list(data.keys())
    
    print(files)
    count=len(files)-1
    count_for_individualfiles=len(files)-1
    print(count)
    %matplotlib inline
    df2 = pd.DataFrame(data)
    print(df2)
    
    if not os.path.isfile('myfile.txt')==True: #if the file does not exist
        df2.to_csv('myfile.txt',index=None, sep='\t', mode='a', header=['nstep',  'time', '        temp', 'press', 'etot', '        ektot', '        eptot', '        bond', '        angle', '        dihed', '        ub', '        imp', '        cmap', 'nb', '        eel', '        vdwaals', '        eelec', '        ehbond', 'restri.', 'ekcmt', 'virial', 'volume', '        density']) #or 
    else: # else it exists so append without writing the header
            df2.to_csv('myfile.txt',index=None, sep='\t', mode='a', header=False)
    
    #OUTPUTTING GRAPHS
    #current_dir=(pathlib.Path().absolute())

    df2=df2.astype(float)
    df2['time']=df2['time'].div(1000)
    for i in range(count):
        print(files[count])
        if files[count]== 'time':
            break
        else:
            #df2['time']=df2['time'].div(1000)
            df2.plot('time',files[count])
            plt.show()
            plt.tight_layout()
            
            #Save figure at 300dpi
            #plt.savefig(str(files[count])+'.jpeg', bbox_inches='tight',dpi=300)
            count=count-1
   
    #ask user if they want to output individual files:
    many_files=input("Do you want to output individual files for each parameter? [y/n]")
    
    if many_files == 'n':
          print('Analyze Files is complete')
    else:
        while count_for_individualfiles>=0: #19
            with open(str(files[count_for_individualfiles])+'.txt', 'w') as f:
                df = pd.DataFrame(data[str(files[count_for_individualfiles])])
                df['time']= data['time']                   #adding time to the files
                df.to_csv(str(files[count_for_individualfiles]) + '.txt', header=[str(files[count_for_individualfiles]),'time'] ,index=None, sep='\t', mode='w') #or r'(directory)
                print('creating file '+str(files[count_for_individualfiles]))
                count_for_individualfiles=count_for_individualfiles-1
        
        print('Analyze Files is complete')
    
    
#Check that programs are installed, if not install them
def install():
    
    required = {'matplotlib', 'pandas'}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        python = sys.executable
        subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)
        print('downlading required packages')
    else:
        print('all required packages are downloaded ')

    
def main():
    install()
    try:
        choice= int(input('How would you like to read your output file(s)?(please indicate by number):'
                  '\n1)To give the direct name of the output file(s)--enter "1" '
                  '\n2)To give the ending extension of your output file(s)[ex).out]--enter "2"\n'))
        #print(type(choice))
        if choice==2:
            ChoiceofOutput(choice)
        else:
            #print('nope')
            GetFiles(choice)
    except TypeError or ValueError:
        print('please indicate by number')
    

if __name__ == '__main__':
    main()
