import time
from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager
from RegexParser import *
import sys
import pandas as pd
from bs4 import BeautifulSoup
import re as re
import numpy as np
from tqdm import tqdm

#Initialize browser
my_parser = RegexParser() 
f=open("tt.txt","w")
driver = webdriver.Chrome(ChromeDriverManager().install())



def main(parameter,dic):


	#Open the URL
	driver.get(parameter)
	time.sleep(2) # Let the user actually see something!


	content = driver.page_source

	soup = BeautifulSoup(content,'html.parser') # getting the html content of this url

	
	gene = parameter.split('/')[3]
	tissue = parameter.split('/')[5].split('#')[0].replace('+',' ')

	TPM_list = []
	for i in soup.findAll('table',attrs={'id':'rnaseq'}):
	
		for j in i.findAll('td', attrs={'class':'center white'}):

			TPM=j.find('b').text
			TPM_list.append(TPM)
			
				
			
	dic2={}
   

	dic[tissue] = TPM_list
	dic2[gene] = dic
	print(dic2)

	name = gene+tissue
	
	np.save(f'/Users/xinyutang/Documents/HPA2/Data/HPA-V8/{name}.npy',dic2)

	
			




if __name__== "__main__":

	url = np.load('/Users/xinyutang/Documents/HPA2/autophagyV8Url_added.npy',allow_pickle=True)
 
	for i in tqdm(url):
		print(i)

		dic={}
		try:
			main(i,dic)
		except Exception:
			print("Oops!  That was no web")
		
		

	
