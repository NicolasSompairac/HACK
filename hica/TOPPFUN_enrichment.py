import requests
import json
import pandas as pd
from json2html import *


def Transform_gene_list_to_dict(gene_list):
    """
    Transform a list of gene ids into a dictionary to serve as API input.
    
    Inputs
    ------
    gene_list: list of strings
        List of strings containing genes ids.
    
    Outputs
    -------
    genes_dict: dict
        Dictionary with key="Symbols" and value=list of genes.
    
    """
    
    genes_dict = {}
    genes_dict["Symbols"] = []
    for gene in gene_list:
        genes_dict["Symbols"].append(gene)
    return genes_dict


def Call_ID_conversion_API(gene_dict):
    """
    Call TOPPGENE API to convert gene ids: https://toppgene.cchmc.org/API/lookup.
    Genes as input can use the following ids: HGNC, ENSEMBL, ENTREZ, UNIPROT, REFSEQ.
    Output ids will be: OfficialSymbol, Entrez, Submitted.
    
    Inputs
    ------
    gene_dict: dict
        Dictionary with key="Symbols" and value=list of genes.
    
    Outputs
    -------
    response: requests.models.Response
        Response from "https://toppgene.cchmc.org/API/lookup" API call.
    
    """
    
    url = "https://toppgene.cchmc.org/API/lookup"
    headers = {'Content-Type': 'text/json'}
    response = requests.post(url,headers=headers,data=json.dumps(gene_dict))
    if response.status_code == 200:
        print("ID conversion success!")
    else:
        print("Something went wrong during conversion... Status code:", response.status_code)
    return response


def Extract_Entrez_id(api_response):
    """
    Extract Entrez IDs from the TOPPGENE lookup API call and store it into
    a dictionary to use as input for TOPPFUN enrichment API.
    
    Inputs
    ------
    api_response: requests.models.Response
        Response from "https://toppgene.cchmc.org/API/lookup" API call.
    
    Outputs
    -------
    entrez_dict: dict
        Dictionary with key="Genes" and value=list of genes in Entrez IDs.
    
    """
    
    entrez_dict = {}
    entrez_dict["Genes"] = []
    for item in api_response.json()["Genes"]:
        entrez_dict["Genes"].append(item["Entrez"])
    return entrez_dict


def Call_enrichment_API(entrez_dict, type_list, Pval=0.05, mingene=1, maxgene=500, maxres=10, correct="FDR"):
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
    entrez_dict: dict
        Dictionary with key="Genes" and value=list of genes in Entrez IDs.
    
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
    data_all["Genes"] = entrez_dict["Genes"]
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
