import requests
import json
import pandas as pd
from json2html import *
import mygene


def Call_ID_conversion_API(genes_list, input_type, output_type, remove_duplicates=True):
	"""
	Call mygene "quarymany()" function to convert gene ids: https://docs.mygene.info/projects/mygene-py/en/latest/#mygene.MyGeneInfo.querymany
	Genes as input can use the following ids: HGNC, ENSEMBL, ENTREZ, UNIPROT, REFSEQ.
	Output ids will be: OfficialSymbol, Entrez, Submitted.
	
	Inputs
	------
	genes_list: list of strings
		List containing gene IDs.

	input_type: string
		Type of input gene IDs.
		For available types, see: https://docs.mygene.info/en/latest/doc/query_service.html#available_fields

	output_type: string
		Type of output gene IDs.
		For available types, see: https://docs.mygene.info/en/latest/doc/data.html#available-fields

	remove_duplicates: Boolean
		If True (default), only first entry IDs (highest matching scores) will be kept.
	
	Outputs
	-------
	response: pandas Dataframe
		Dataframe containing query IDs as index.
	
	"""
	
	#Initiate mygene
	mg = mygene.MyGeneInfo()
	#Perform Conversion
	conversion_df = mg.querymany(genes_list, scopes=input_type, fields=output_type, verbose=False, as_dataframe=True)
	print("ID conversion success!")

	if remove_duplicates:
			return conversion_df[~conversion_df.index.duplicated(keep='first')]
	else:
		return conversion_df

def Extract_Entrez_id(api_response_df):
	"""
	Extract Entrez IDs from the mygene conversion API call and store it into
	a list to use as input for TOPPFUN enrichment API.
	
	Inputs
	------
	api_response_df: pandas Dataframe
		Dataframe from mygene "querymany()" function containing Entrez IDs.
	
	Outputs
	-------
	entrez_list: list
		List of genes in Entrez IDs.
	
	"""
	
	entrez_list = list(map(int, list(api_response_df["entrezgene"])))
	
	return entrez_list


def Call_enrichment_API(entrez_list, type_list=None, Pval=0.05, mingene=1, maxgene=500, maxres=10, correct="FDR"):
	"""
	Call TOPPGENE API to detect functional enrichments of a gene list: https://toppgene.cchmc.org/API/enrich.
	Enrichment can be performed on the following features: GeneOntologyMolecularFunction
														   GeneOntologyBiologicalProcess
														   GeneOntologyCellularComponent
														   HumanPheno
														   MousePheno
														   Domain
														   Pathway
														   Pubmed
														   Interaction
														   Cytoband
														   TFBS
														   GeneFamily
														   Coexpression
														   CoexpressionAtlas
														   GeneFamily
														   Computational
														   MicroRNA
														   Drug
														   Disease
	
	Inputs
	------
	entrez_list: list
		List of genes in Entrez IDs.
	
	type_list: list of strings
		List of features to perform enrichment tests.
	
	Pval: float
		P-value maximal threshold to accept enrichments.
	
	mingene: int
		Minimal number of gene hits to accept in feature categories.
	
	maxgene: int
		Maximal number of gene hits to accept in feature categories.
	
	maxres: int
		Number of top enrichments to show for each feature.
	
	correct: str {'none', 'FDR', 'Bonferroni'}
		P-value correction methods. The FDR method is the Benjamini and Hochberg method.
	
	Outputs
	-------
	response: requests.models.Response
		Response from "https://toppgene.cchmc.org/API/enrich" API call.
	
	"""

	if type_list is None :
			type_list = ["GeneOntologyMolecularFunction","GeneOntologyBiologicalProcess","GeneOntologyCellularComponent","HumanPheno","MousePheno",
						 "Domain","Pathway","Pubmed","Interaction","Cytoband","TFBS","GeneFamily","Coexpression","CoexpressionAtlas","GeneFamily",
						 "Computational","MicroRNA","Drug","Disease"]
	
	url = "https://toppgene.cchmc.org/API/enrich"
	headers = {'Content-Type': 'text/json'}
	parameters = {}
	parameters["Categories"] = []
	for type_id in type_list:
		parameters["Categories"].append({"Type":type_id,
										"Pvalue":Pval,
										"MinGenes":mingene,
										"MaxGenes":maxgene,
										"MaxResults":maxres,
										"Correction":correct})
	data_all = {}
	data_all["Genes"] = entrez_list
	data_all["Categories"] = parameters["Categories"]
	response = requests.post(url,headers=headers,data=json.dumps(data_all))
	if response.status_code == 200:
		print("Enrichment analysis success!")
	else:
		print("Something went wrong during enrichment... Status code:", response.status_code)
	return response


def Clean_enrichment_genes(api_response):
	"""
	Clean gene IDs from the TOPPGENE enrichment API call reponse for a prettier output.
	
	Inputs
	------
	api_response: requests.models.Response
		Response from "https://toppgene.cchmc.org/API/enrich" API call.
	
	Outputs
	-------
	clean_json: list of dict
		List containing dictionaries, simulating the TOPPGENE enrichment API response,
		containing only hit Gene Symbols.
	
	"""
	
	clean_json = []
	#Clean_json["Annotations"] = []
	for element in api_response.json()["Annotations"]:
		#print(element)
		gene_symbol_list = []
		for gene in element["Genes"]:
			#print(gene, gene["Symbol"])
			gene_symbol_list.append(gene["Symbol"])
		#element["Gene_Symbol"] = gene_symbol_list
		element["Gene_Symbol"] = ','.join(gene_symbol_list)
		element.pop("Genes", None)
		#print(element)
		clean_json.append(element)

	return clean_json


def Export_enrichment_results(cleaned_json, filetype, filename):
	"""
	Export TOPPFUN enrichment results in a HTML or Excel file.
	
	Inputs
	------
	cleaned_json: list
		Cleaned list of the response from "https://toppgene.cchmc.org/API/enrich" API call.
	
	filetype: str {'html', 'xslx'}
		Output file type.
	
	filename: str
		Root name of the output file, without the extension.
	
	Outputs
	-------
	None
	
	"""
	
	if filetype == "html":
		html_format = json2html.convert(json = cleaned_json)
		with open(filename+".html", "w") as outfile:
			outfile.write(html_format)
		print("Exported enrichment result in HTML!")
	elif filetype == "xlsx":
		df = pd.DataFrame(cleaned_json)
		df.to_excel(filename+".xlsx")
		print("Exported enrichment result in Excel!")
	else:
		"Wrong type input. Please use 'html' or 'xlsx'."
	
	return



# Example list of genes related to immune receptor activity.
Gene_query_list = ["HLA-DOB","HLA-DQA1","HLA-DRA","FCER1G","PIGR","CCR4","CCR6","CCR7","IL27RA","FPR3","CR1","KLRK1","CSF2RB","IL7R","CXCR5","IL22RA2","CXCR4","C3AR1","CCR10"]
# Parameter for ID conversion
Input_type = 'symbol'
Output_type = 'entrezgene'
# Convert gene ids to Entrez
Conversion_API_response = Call_ID_conversion_API(Gene_query_list,Input_type,Output_type)
# Extract Entrez ids from the conversion API response
Entrez_ids_list = Extract_Entrez_id(Conversion_API_response)
# Choose the list of features to perform enrichment against
Enrichment_types = ["GeneOntologyMolecularFunction",
					"GeneOntologyBiologicalProcess",
					"GeneOntologyCellularComponent",
					"HumanPheno",
					"MousePheno",
					"Domain",
					"Pathway",
					"Pubmed",
					"Interaction",
					"Cytoband",
					"TFBS",
					"GeneFamily",
					"Coexpression",
					"CoexpressionAtlas",
					"GeneFamily",
					"Computational",
					"MicroRNA",
					"Drug",
					"Disease"]
# Enrichment parameters
Pval = 0.05
Min_genes = 1
Max_genes = 500
Max_results = 5
Correction = "FDR"
# Call TOPPFUN enrichment API
Enrichment_API_response = Call_enrichment_API(Entrez_ids_list,Enrichment_types,
											  Pval,Min_genes,Max_genes,Max_results,
											  Correction)
# Clean the enrichment response JSON structure for pretier exports
Enrichment_API_response_cleaned = Clean_enrichment_genes(Enrichment_API_response)
# Export enrichment analysis into a file
File_type = "html" # use html or xlsx
File_name = "Enrichment_example"
Export_enrichment_results(Enrichment_API_response_cleaned, File_type, File_name)
File_type = "xlsx" # use html or xlsx
Export_enrichment_results(Enrichment_API_response_cleaned, File_type, File_name)
