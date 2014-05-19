#!/usr/bin/env python

import re, json, os

def removeComments(string): 	#Removes comments denoted by '//'
    string = re.sub(re.compile("//.*" ) ,"" ,string)
    return string

def ImportJSON(filename):
	jsonin = open(filename,"r")
	jsontemp = open("inputs_temp.json","w")

	for line in jsonin:
		jsontemp.write(removeComments(line))	#Parse through the json file, removes problematic comments
	jsontemp.write("\n")
	jsonin.close()
	jsontemp.close()

	json_data = open("inputs_temp.json","r")	#Open the copy of the file sans comments
	params = json.load(json_data)		#Load the json file into a dictionary
	json_data.close()
	os.system("rm inputs_temp.json")
	return params
