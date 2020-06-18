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

	
	gene = parameter.split('/')[3].split('-')[1]
	
	gene_dic={}
	celline_dic={}
	block_matches = my_parser.parse("<g class=\"bar_g\"(.*?)<br>Category",driver.page_source)
	count=0
	l=[]



	for match in block_matches:
		count+=1
	
		Celline= my_parser.parse('title=\"<b>(.*?)<\/b>',match)[0]
		Organ =my_parser.parse('Organ: (.*?)<br>',match)[0]

		celline_dic[Celline] =Organ

		print(count)

	
	gene_dic[gene]=celline_dic
	print(gene_dic)
	
	
		
	with open(f'/Users/xinyutang/Documents/HPA2/Data/HPA-V9/{gene}.json', 'w') as fp:
		json.dump(gene_dic, fp)
			




if __name__== "__main__":

	url = np.load('/Users/xinyutang/Documents/HPA2/crawler/cellurl.npy',allow_pickle=True)
 
	for i in tqdm(url):
		print(i)

		dic={}
		try:
			main(i,dic)
		except Exception:
			print("Oops!  That was no web")
		
		

	

	
