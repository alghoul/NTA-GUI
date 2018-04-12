# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 10:11:13 2017

@author: HALGhoul
"""

import pandas as pd
import numpy as np
import re
import os
from pandasql import sqldf
from operator import itemgetter
from itertools import groupby
from difflib import get_close_matches
from difflib import SequenceMatcher
from itertools import izip
#REP_NUM = 3
HBR = 3.0 # High_Blank_Ratio condition
HMR = 1.5 # High_Mid_Ratio condition
SCORE = 90 # formula match is 90


def fix_names(df,index): # parse the Dataframe into a numpy array
    	#df.columns = df.columns.str.replace(': Log2','') #log specific code
    	df.columns = df.columns.str.replace(' ','_')
    	df.columns = df.columns.str.replace('\([^)]*\)','')
    	df['Compound'] = df['Compound'].str.replace("\Esi.*$","")
    	#df.drop(['CompositeSpectrum','Compound_Name'],axis=1)
    	df.drop(['Compound_Name'],axis=1)
    	Headers = parse_headers(df,index)
    	Abundance = [item for sublist in Headers for item in sublist if len(sublist)>1]	
	Samples= [x for x in Abundance]
	NewSamples = common_substrings(Samples)
	for i in range(len(Samples)):
		df.rename(columns = {Samples[i]:NewSamples[i]},inplace=True)	
    	#df = df
    	return df





def read_data(file,index):  # read a csv file into a DataFrame
	ext = os.path.splitext(file)[1]
	print ext
	if ext == '.tsv':
	    	df = pd.read_csv(file,sep='\t',comment='#',na_values= 1 | 0)
	if ext == '.csv':
	    	df = pd.read_csv(file,comment='#',na_values= 1 | 0)
    	df = fix_names(df,index)
    	return df





def differences(s1,s2): #find the number of different characters between two strings (headers)
    	s1 = re.sub(re.compile(r'\([^)]*\)'),'',s1)
    	s2 = re.sub(re.compile(r'\([^)]*\)'),'',s2)    
    	count = sum(1 for a, b in zip(s1, s2) if a != b) #+ abs(len(s1) - len(s2))
    	return count




def formulas(df):
    	formulas = df.loc[df['For_Dashboard_Search'] == '1','Compound'].values #only features flagged for Dashboard search
    	print formulas
    	return formulas




def parse_headers(df,index): #group headers into a group of samples
    	#global df
    	headers = [[],[]]
    	headers[index] = df.columns.values.tolist()
    	countS=0
    	countD=0
    	new_headers = [[],[]]
    	New_Headers = [None,None]
    	Headers = [None,None]
    	groups = [None,None]
    	for s in range(0,len(headers[index])-1):
        	#print headers[s],headers[s+1],list(set(str(headers[s])) - set(str(headers[s+1])))
        	if differences(str(headers[index][s]),str(headers[index][s+1])) < 2:
            		countS += 1
            		#countD = 0
            		#print "likely the same sample " 
        	if differences(str(headers[index][s]),str(headers[index][s+1])) >= 2:
            		countD += 1
            		countS = countS + 1
            		#print "These are different "
		if "_Flags" in headers[index][s]:
	    		break
        	new_headers[index].append([headers[index][countS],countD])
        	new_headers[index].sort(key = itemgetter(1))
        	groups[index] = groupby(new_headers[index], itemgetter(1))
        	New_Headers[index] = [[item[0] for item in data] for (key, data) in groups[index]] 
    	Headers[index] = New_Headers[index]
    	return Headers[index]





def score(df): # Get that Sneaky score from Annotations.
    	regex = "^.*=(.*) \].*$" # a regex to find the score looking for a pattern of "=something_to_find ]" 
	if "Annotations" in df:
	    	df['Score'] = df.Annotations.str.extract(regex,expand=True).astype('float64')
    	return df





def statistics(df,index): # calculate Mean,Median,STD,CV for every feature in a sample of multiple replicates
	Abundance = [[],[]]
	NewAbundance = [[],[]]
	Blanks = [[],[]]
	NewBlanks = [[],[]]
	Samples = [[],[]]
	NewSamples = [[],[]]	
    	Headers = [0,0]
    	NewHeaders = [0,0]
    	headers = [0,0]
    	Headers[index] = parse_headers(df,index)
    	Abundance[index] = [item for sublist in Headers[index] for item in sublist if len(sublist)>1]
    	df  = score(df) 

	# Do some statistical acrobatics
    	headers[index] = ['Compound','Ionization_Mode','Score','Mass','Retention_Time','Frequency'] + Abundance[index]
    	df = df[headers[index]]
	print Headers[index] #stopped here before my optometrist appointment
    	for list in Headers[index]:
            	REP_NUM = len(list)
            	if REP_NUM > 1:
                	for i in range(0,REP_NUM):
                    	# the match part finds the indices of the largest common subtring between two strings
                    		match = SequenceMatcher(None, list[i], list[i+1]).find_longest_match(0, len(list[i]),0, len(list[i+1]))
                    		df['Mean_'+ str(list[i])[match.a:match.a +  match.size]] = df[list[i:i + REP_NUM]].mean(axis=1).round(0)
                    		df['Median_'+ str(list[i])[match.a:match.a +  match.size]] = df[list[i:i + REP_NUM]].median(axis=1,skipna=True).round(0) 
                    		df['STD_'+ str(list[i])[match.a:match.a +  match.size]] = df[list[i:i + REP_NUM]].std(axis=1,skipna=True).round(0)
                    		df['CV_'+ str(list[i])[match.a:match.a +  match.size]] = (df['STD_'+ str(list[i])[match.a:match.a +  match.size]]/df['Mean_'+ str(list[i])[match.a:match.a +  match.size]]).round(2)            
                    		df['N_Abun_'+ str(list[i])[match.a:match.a +  match.size]] = df[list[i:i + REP_NUM]].count(axis=1).round(0)
                    		#print list[i][match.a:match.a +  match.size]
                    		break
    	df.sort_values(['Mass','Retention_Time'],ascending=[True,True],inplace=True)    
    	#df.to_csv('input-updated.csv', index=False)
    	return df



def Blank_Subtract(df,index):
	Abundance = [[],[]]
    	Headers = [0,0]
	Blanks = [[],[]]
    	Headers[index] = parse_headers(df,index)
    	Abundance[index] = [item for sublist in Headers[index] for item in sublist if len(sublist)>1]

	# On with the agony of subtracting the MB median from Samples
	Blanks[index] = df.columns[df.columns.str.contains(pat ='MB_|blank|blanks')].tolist()

	df['Median_ALLMB'] = df[Blanks[index]].median(axis=1,skipna=True).round(0).fillna(0)
	df[Abundance[index]] = df[Abundance[index]].sub(df['Median_ALLMB'],axis=0) #subtract the median of MBs from every Sample median
	df[Abundance[index]] = df[Abundance[index]].clip(lower=0).replace({0:np.nan})
    	return df


 
def check_feature_tracers(df,index,Mass_Difference,Retention_Difference): #a method to query and save the features with tracers criteria
    	df_sql = [None,None]
    	Statistics = [[],[]]
    	b_Statistics = [[],[]]
    	q = [None,None]
    	df1 = df    
    	df2 = pd.read_csv("Tracers_Table_ for_SRM2585_20170524.csv",comment='#',na_values= 1 | 0)
    	#read_data("Tracers_Table_ for_SRM2585_20170524.csv")
    	Statistics[index] = df.columns[df.columns.str.contains(pat ='N_|CV_|Mean_|Median_|STD_')].tolist()
    	print df2
    	b_Statistics[index] = ["b." + B for B in Statistics[index]]
    	q[index] ="""select	a.*, b.Mass as Observed_Mass,
				 b.Retention_Time as Observed_Retention_Time,""" + " , ".join([b + " as " + a  for b,a in zip(b_Statistics[index],Statistics[index])]) + """ from df2 as a left join df1 as b where abs((a.Monoisotopic_Mass-b.Mass)/a.Monoisotopic_Mass)*(1000000)<=""" + str(Mass_Difference) + """ and abs(a.Retention_Time-b.Retention_Time)<=""" + str(Retention_Difference) + """ and a.Ionization_Mode = b.Ionization_Mode;"""
    	df_sql[index] = sqldf(q[index],locals())         
    	#df_sql.to_csv('input_after_tracers.csv', index=False)
    	return df_sql[index]





def clean_features(df,index): # a method that drops rows based on conditions
    	Abundance=[[],[]]
    	Abundance[index] =  df.columns[df.columns.str.contains(pat ='N_Abun_')].tolist()   
    	#for header in Abundance:
		#    df = df.drop(df[df[header] < 2].index) # drop rows with n_abundance_high <2
    	'''
    	CV=list()
    	CV =  df.columns[df.columns.str.contains(pat ='CV_')].tolist()  
    	for cv in CV: 
	    	df = df.drop(df[df[cv] != df[cv]].index) '''    
    	Median=[[],[]]
    	Median_MB = [[],[]]
    	N_Abun_MB = [[],[]]
    	Median[index] =  df.columns[df.columns.str.contains(pat ='Median_')].tolist() 
    	Median_MB[index] = [md for md in Median[index] if 'MB' in md]
    	N_Abun_MB[index] = [N for N in Abundance[index] if 'MB' in N]
    	for median in Median[index]:
	    	df = df.drop(df[(df[median]/df[Median_MB[index][0]] < HBR) | (df[N_Abun_MB[index][0]] != 0)].index)   
    	#df = df.drop(df[df.Score != df.Score].index)
    	return df
  



    
def flags(df): # a method to develop required flags
    	df['Neg_Mass_Defect'] = np.where((df.Mass - df.Mass.round(0)) < 0 , '1','0')
    	df['Halogen'] = np.where(df.Compound.str.contains('F|l|r|I'),'1','0')
    	df['Formula_Match'] = np.where(df.Score != df.Score,'0','1') #check if it does not have a score
    	df['Formula_Match_Above90'] = np.where(df.Score >= SCORE,'1','0')
    	df['X_NegMassDef_Below90'] = np.where(((df.Score < SCORE) & (df.Neg_Mass_Defect == '1') & (df.Halogen == '1')),'1','0')
    	df['For_Dashboard_Search'] = np.where(((df.Formula_Match_Above90 == '1') | (df.X_NegMassDef_Below90 == '1')) , '1', '0') 
    	df.sort_values(['Formula_Match','For_Dashboard_Search','Formula_Match_Above90','X_NegMassDef_Below90'],ascending=[False,False,False,False],inplace=True) 
    	#df.to_csv('input-afterflag.csv', index=False) 
    	#print df1 
    	df.sort_values('Compound',ascending=True,inplace=True)
    	return df



def match_headers(list1=None,list2=None):
	match = None
	string_match = list()
	print len(list1), len(list2)
	for i in range(len(list2)):
		print list1[i] + "  ,  " + list2[i] 
		string_match.append("".join([list2[i][j] for j, (a,b) in enumerate(izip(list1[i],list2[i])) if a == b]))	
	print len(string_match)
	return string_match	
	

def common_substrings(ls=[]):
	match  = SequenceMatcher(None,ls[0],ls[len(ls)-1]).find_longest_match(0,len(ls[0]),0,len(ls[len(ls)-1]))
	common = ls[0][match.a: match.a + match.size]
	print " ********* " + common
	lsnew = list()
	for i in range(len(ls)):
		lsnew.append(ls[i].replace(common,''))
	#print ls
	return lsnew
	


def combine(df1,df2):
	#Headers = [[],[]]
    	#Headers[0] = parse_headers(df1,0)
    	#Headers[1] = parse_headers(df2,1)
	print "##############"
	Abundance=[[],[]]
	Abundance[0] = df1.columns.values.tolist()
	Abundance[1] = df2.columns.values.tolist()	 	
	new_headers = match_headers(Abundance[0],Abundance[1])
	#print len(df1.columns.values.tolist())
	for i in range(len(Abundance[0])):
		#print (Abundance[0][i],Abundance[1][i])
		df1.rename(columns = {Abundance[0][i]:new_headers[i]},inplace=True)
		df2.rename(columns = {Abundance[1][i]:new_headers[i]},inplace=True)
	#print df1.columns.values.tolist()
	print " ||||___|||| - - - - - - "
	#print df2.columns.values.tolist()	
	dfc = pd.concat([df1,df2]).reset_index()
	columns = dfc.columns.values.tolist()
	dfc = pd.merge(dfc,df2,suffixes=['','_x'],on='Compound',how='left')
	dfc = pd.merge(dfc,df1,suffixes=['','_y'],on='Compound',how='left')

	# create new flags
	dfc = dfc.drop_duplicates(subset=['Compound','Mass','Retention_Time','Score'])
	dfc['Both_Modes'] = np.where(((abs(dfc.Mass_x-dfc.Mass_y)<=0.005) & (abs(dfc.Retention_Time_x-dfc.Retention_Time_y)<=1)),'1','0')
	dfc['N_Compound_Hits'] = dfc.groupby('Compound')['Compound'].transform('size')
	Median_list =  dfc.columns[dfc.columns.str.contains(pat ='Median_')].tolist()
	#print Median_list 	
	dfc['N_Abun_Samples'] = dfc[Median_list].count(axis=1)
	dfc['Median_Abun_Samples'] = dfc[Median_list].median(axis=1,skipna=True).round(0)
	dfc['One_Mode_No_Isomers'] = np.where(((dfc.Both_Modes == '0') & (dfc.N_Compound_Hits == 1)),'1','0')
	dfc['One_Mode_Isomers'] = np.where(((dfc.Both_Modes == '0') & (dfc.N_Compound_Hits > 1)),'1','0')
	dfc['Two_Modes_No_Isomers'] = np.where(((dfc.Both_Modes == '1') & (dfc.N_Compound_Hits == 2)),'1','0')
	dfc['Two_Modes_Isomers'] = np.where(((dfc.Both_Modes == '1') & (dfc.N_Compound_Hits > 2)),'1','0')
	dfc['Est_Chem_Count'] = None #Default to non-type
	dfc['Est_Chem_Count'][dfc['One_Mode_No_Isomers'] == '1'] = 1
	dfc['Est_Chem_Count'][dfc['One_Mode_Isomers'] == '1'] = dfc['N_Compound_Hits']
	dfc['Est_Chem_Count'][(dfc['Two_Modes_No_Isomers'] == '1') | (dfc['Two_Modes_Isomers'] == '1')] = dfc['N_Compound_Hits']/2	
	columns.extend(('Both_Modes','N_Compound_Hits','N_Abun_Samples','Median_Abun_Samples','One_Mode_No_Isomers','One_Mode_Isomers','Two_Modes_No_Isomers',
			'Two_Modes_Isomers','Est_Chem_Count'))
	dfc = dfc[columns].sort_values(['Compound'],ascending=[True])

	#dft.reset_index() 
	#dft.dropna(inplace=True)
	return dfc
	



def reduce(df,index):
    	Abundance = [[],[]]
    	Headers = [0,0]
    	Headers[index] = parse_headers(df,index)
    	Abundance[index] = [item for sublist in Headers[index] for item in sublist if len(sublist)>2] 	    
    	df.drop(Abundance[index],axis=1,inplace=True)
    	return df





def duplicates(df,index):
    	df1 = df
    	Abundance = [[],[]]
    	q = [None,None]
    	Abundance[index] = [item for sublist in parse_headers(df,index) for item in sublist if len(sublist)>1]
    	#print Abundance[index]
    	#df1[Abundance[index]] = df1[Abundance[index]].astype('float64').apply(np.exp2).round(0).replace(1,np.NaN) #log specific code
    	q[index]="""select a.*, b.Mass as match_Mass, b.Retention_Time as match_Retention_Time
		from df1 as a, df1 as b where (abs(a.Mass-b.Mass)<=0.005 
		and abs(a.Retention_Time-b.Retention_Time)<=0.05 and (a.Compound == b.Compound or a.Compound LIKE '%@%')) order by a.Mass, b.Mass; """
    	df2 = sqldf(q[index],locals()) 
    	df2.sort_values(['Mass','Retention_Time'],ascending=[True,True],inplace=True)
    	df2['Mass_Rounded'] = df2['Mass'].round(2)
    	df2['Retention_Time_Rounded'] = df2['Retention_Time'].round(1)
    	df2.to_csv('Testing-duplicates-Algorithm-before.csv', index=False)
    	for i in Abundance:
    	    	df2[i] = df2[i].fillna(df2.groupby(['Mass_Rounded','Retention_Time_Rounded'])[i].transform(max))
    	df2.to_csv('Testing-duplicates-Algorithm-after.csv', index=False)
    	df2 = df2.drop_duplicates(subset=['Mass_Rounded','Retention_Time_Rounded'], keep="last")
    	df = df2
    	return df



def MPP_Ready(dft, directory='',file=''):
	dft = dft.rename(columns = {'Compound':'Formula','Retention_Time':'RT'})
	dft['Compound Name'] = dft['Formula']
	dft['CAS ID'] = ""
    	Headers = parse_headers(dft,0)
    	Abundance = [item for sublist in Headers for item in sublist if len(sublist)>2]
	Blanks = dft.columns[dft.columns.str.contains(pat ='MB_')].tolist()
	Samples = [x for x in Abundance if x not in Blanks]
	NewSamples = common_substrings(Samples)
	for i in range(len(Samples)):
		dft.rename(columns = {Samples[i]:NewSamples[i]},inplace=True)
	#columns = dft.columns.values.tolist()
	#dft = dft.reindex(columns=Columns)
	#print dft
	#dft.to_csv(directory+'/'+file+'_MPP_Ready.csv', index=False)
	dft = dft[['Formula','Compound Name','CAS ID','Mass','RT'] + NewSamples]	
	dft.to_csv(directory+'/'+file+'_MPP_Ready.csv', index=False)	
	return dft







#check_feature_tracers(read_data("House_Dust_Negative_MPP_output.csv"),read_data("Tracers_Table_ for_SRM2585_20170524.csv"))


