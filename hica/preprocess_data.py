import numpy as np
import scanpy as sc
import pandas as pd
import os

def Load_data(foldername, filename, delim="\t", cache_state=True):
	
	df = pd.read_csv(foldername+"/"+filename,sep=delim,header=0,index_col=0).transpose()
	adata = sc.AnnData(df)
	
	return adata

def Preprocess_data(adata, To_Log=True, Keep_Variable_Genes=0, Min_cells_prop=0, Center=True, Remove_Duplicates=True):
	
	if To_Log:
		sc.pp.log1p(adata)
	if Keep_Variable_Genes:
		var = np.var(adata.X,axis=0)
		inds = np.flip(np.argsort(var))
		ind_genes = inds[0:Keep_Variable_Genes]
		if 0 in var[ind_genes]:
			ind_first_zero = np.argwhere(var[ind_genes]==0)[0][0]
			ind_genes = ind_genes[0:ind_first_zero]
		adata = adata[:,ind_genes]
	if Min_cells_prop:
		sc.pp.filter_genes(adata, min_cells = int((Min_cells_prop*len(adata.X))), inplace = True)
	df = adata.copy().to_df()
	if Remove_Duplicates:
		df = df.loc[:,~df.columns.duplicated()]
	if Center:
		df = df.apply(lambda x: x-x.mean())
	
	return df

def Save_data(df, foldername, job):
	
	df_T = df.transpose()
	df_T.to_csv(foldername+"/"+job+"_proc.txt", sep='\t')

	# Extract genes in a separate file
	with open(foldername+"/"+job+"_genes.txt", 'w') as outfile:
		outfile.write("Genes\n")
		for gene in df_T.index.values:
			outfile.write(str(gene)+"\n")
	
	return

def Load_as_numpy(foldername, job):
	
	full_path = foldername+"/"+job
	SC_data = sc.read(full_path+"_proc.txt",
					cache = True, first_column_names = True, delimiter = "\t").transpose()
	np.save(file=full_path+"_proc", arr=SC_data.X.transpose())
	SC_data = np.load(full_path+'_proc.npy')
	
	return SC_data
