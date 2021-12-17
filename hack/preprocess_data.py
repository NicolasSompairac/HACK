import numpy as np
import scanpy as sc
import pandas as pd
import os

def Preprocess_data_new(foldername, filename, job, delim="\t", To_Log=True, Keep_Variable_Genes=0, Min_cells_prop=0, Center=True, Remove_Duplicates=True):

	adata = sc.read_text(filename=foldername+"/"+filename, delimiter=delim,first_column_names=True).transpose()

	# Preprocessing
	print("Performing preprocessing")
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
	df = adata.to_df()
	if Remove_Duplicates:
		df = df.loc[:,~df.columns.duplicated()]
	if Center:
		df = df.apply(lambda x: x-x.mean())
	

	full_path = foldername+"/"+job
	# Extract genes in a separate file
	with open(full_path+"_genes.txt", 'w') as outfile:
		outfile.write("Genes\n")
		for gene in df.columns.values:
			outfile.write(str(gene)+"\n")

	# Extract samples in a separate file
	with open(full_path+"_samples.txt", 'w') as outfile:
		outfile.write("Samples\n")
		for sample in df.index.values:
			outfile.write(str(sample)+"\n")

	# Save processed data in numpy
	print("Saving result as numpy file")
	adata = sc.AnnData(df)
	np.save(file=full_path+"_proc", arr=adata.X.transpose())
	print("Done!")

	return adata.X.transpose()

