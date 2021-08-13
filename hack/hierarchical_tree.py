import pandas as pd
import numpy as np
import scanpy as sc
import scipy.stats
import networkx as nx
import os

def Create_full_graph(datafolder, biton_folder, job, Minimum_decomposition, Maximum_decomposition):
	
	# Load data genes
	data = sc.read(datafolder+"/"+job+"_genes.txt",
							cache = False, first_column_names = True).transpose()

	print('Preparing metagenes...')
	# Load all Biton metagenes
	Biton_M2_GC_CONTENT = pd.read_csv(biton_folder+'/Biton_M2_GC_CONTENT.txt', sep='\t')
	Biton_M3_SMOOTH_MUSCLE = pd.read_csv(biton_folder+'/Biton_M3_SMOOTH_MUSCLE.txt', sep='\t')
	Biton_M4_MITOCHONRIA_TRANSLATION = pd.read_csv(biton_folder+'/Biton_M4_MITOCHONRIA_TRANSLATION.txt', sep='\t')
	Biton_M5_INTERFERON = pd.read_csv(biton_folder+'/Biton_M5_INTERFERON.txt', sep='\t')
	Biton_M6_BASALLIKE = pd.read_csv(biton_folder+'/Biton_M6_BASALLIKE.txt', sep='\t')
	Biton_M7_CELLCYCLE = pd.read_csv(biton_folder+'/Biton_M7_CELLCYCLE.txt', sep='\t')
	Biton_M8_IMMUNE = pd.read_csv(biton_folder+'/Biton_M8_IMMUNE.txt', sep='\t')
	Biton_M9_UROTHELIALDIFF = pd.read_csv(biton_folder+'/Biton_M9_UROTHELIALDIFF.txt', sep='\t')
	Biton_M12_MYOFIBROBLASTS = pd.read_csv(biton_folder+'/Biton_M12_MYOFIBROBLASTS.txt', sep='\t')
	Biton_M13_BLCAPATHWAYS = pd.read_csv(biton_folder+'/Biton_M13_BLCAPATHWAYS.txt', sep='\t')
	Biton_M14_STRESS = pd.read_csv(biton_folder+'/Biton_M14_STRESS.txt', sep='\t')

	# Find correspondence between data genes and Biton genes
	#Biton_M2_GC_CONTENT
	sign_M2_idx = Biton_M2_GC_CONTENT.index.get_indexer(Biton_M2_GC_CONTENT.index.intersection(data.var.index,sort=None))
	data_M2_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M2_GC_CONTENT.index,sort=None))
	#Biton_M3_SMOOTH_MUSCLE
	sign_M3_idx = Biton_M3_SMOOTH_MUSCLE.index.get_indexer(Biton_M3_SMOOTH_MUSCLE.index.intersection(data.var.index,sort=None))
	data_M3_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M3_SMOOTH_MUSCLE.index,sort=None))
	#Biton_M4_MITOCHONRIA_TRANSLATION
	sign_M4_idx = Biton_M4_MITOCHONRIA_TRANSLATION.index.get_indexer(Biton_M4_MITOCHONRIA_TRANSLATION.index.intersection(data.var.index,sort=None))
	data_M4_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M4_MITOCHONRIA_TRANSLATION.index,sort=None))
	#Biton_M5_INTERFERON
	sign_M5_idx = Biton_M5_INTERFERON.index.get_indexer(Biton_M5_INTERFERON.index.intersection(data.var.index,sort=None))
	data_M5_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M5_INTERFERON.index,sort=None))
	#Biton_M6_BASALLIKE
	sign_M6_idx = Biton_M6_BASALLIKE.index.get_indexer(Biton_M6_BASALLIKE.index.intersection(data.var.index,sort=None))
	data_M6_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M6_BASALLIKE.index,sort=None))
	#Biton_M7_CELLCYCLE
	sign_M7_idx = Biton_M7_CELLCYCLE.index.get_indexer(Biton_M7_CELLCYCLE.index.intersection(data.var.index,sort=None))
	data_M7_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M7_CELLCYCLE.index,sort=None))
	#Biton_M8_IMMUNE
	sign_M8_idx = Biton_M8_IMMUNE.index.get_indexer(Biton_M8_IMMUNE.index.intersection(data.var.index,sort=None))
	data_M8_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M8_IMMUNE.index,sort=None))
	#Biton_M9_UROTHELIALDIFF
	sign_M9_idx = Biton_M9_UROTHELIALDIFF.index.get_indexer(Biton_M9_UROTHELIALDIFF.index.intersection(data.var.index,sort=None))
	data_M9_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M9_UROTHELIALDIFF.index,sort=None))
	#Biton_M12_MYOFIBROBLASTS
	sign_M12_idx = Biton_M12_MYOFIBROBLASTS.index.get_indexer(Biton_M12_MYOFIBROBLASTS.index.intersection(data.var.index,sort=None))
	data_M12_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M12_MYOFIBROBLASTS.index,sort=None))
	#Biton_M13_BLCAPATHWAYS
	sign_M13_idx = Biton_M13_BLCAPATHWAYS.index.get_indexer(Biton_M13_BLCAPATHWAYS.index.intersection(data.var.index,sort=None))
	data_M13_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M13_BLCAPATHWAYS.index,sort=None))
	#Biton_M14_STRESS
	sign_M14_idx = Biton_M14_STRESS.index.get_indexer(Biton_M14_STRESS.index.intersection(data.var.index,sort=None))
	data_M14_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M14_STRESS.index,sort=None))

	print('Calculating correlations...')
	Biton_sign_dict = {}
	for comp in range(Minimum_decomposition,Maximum_decomposition+1):
		X = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(comp)+".npy")
		#print("X", X.shape[1])

		for i in range(X.shape[1]):

			Biton_sign_dict["M2_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M2_idx],Biton_M2_GC_CONTENT["M2_GC_CONTENT"][sign_M2_idx])
			Biton_sign_dict["M3_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M3_idx],Biton_M3_SMOOTH_MUSCLE["M3_SMOOTH_MUSCLE"][sign_M3_idx])
			Biton_sign_dict["M4_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M4_idx],Biton_M4_MITOCHONRIA_TRANSLATION["M4_MITOCHONRIA_TRANSLATION"][sign_M4_idx])
			Biton_sign_dict["M5_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M5_idx],Biton_M5_INTERFERON["M5_INTERFERON"][sign_M5_idx])
			Biton_sign_dict["M6_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M6_idx],Biton_M6_BASALLIKE["M6_BASALLIKE"][sign_M6_idx])
			Biton_sign_dict["M7_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M7_idx],Biton_M7_CELLCYCLE["M7_CELLCYCLE"][sign_M7_idx])
			Biton_sign_dict["M8_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M8_idx],Biton_M8_IMMUNE["M8_IMMUNE"][sign_M8_idx])
			Biton_sign_dict["M9_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M9_idx],Biton_M9_UROTHELIALDIFF["M9_UROTHELIALDIFF"][sign_M9_idx])
			Biton_sign_dict["M12_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M12_idx],Biton_M12_MYOFIBROBLASTS["M12_MYOFIBROBLASTS"][sign_M12_idx])
			Biton_sign_dict["M13_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M13_idx],Biton_M13_BLCAPATHWAYS["M13_BLCAPATHWAYS"][sign_M13_idx])
			Biton_sign_dict["M14_"+str(comp)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][data_M14_idx],Biton_M14_STRESS["M14_STRESS"][sign_M14_idx])

	print("Storing full graph as a file...")
	# Create a folder to store graphs
	if not os.path.exists('Graphs'):
		os.makedirs('Graphs')
	with open('Graphs/'+job+"_full_correlation_graph.txt", 'w') as outfile:

		outfile.write("Source\tTarget\tCorrelation\tBiton_M2_GC_CONTENT\tBiton_M3_SMOOTH_MUSCLE\tBiton_M4_MITOCHONRIA_TRANSLATION\tBiton_M5_INTERFERON\tBiton_M6_BASALLIKE\tBiton_M7_CELLCYCLE\tBiton_M8_IMMUNE\tBiton_M9_UROTHELIALDIFF\tBiton_M12_MYOFIBROBLASTS\tBiton_M13_BLCAPATHWAYS\tBiton_M14_STRESS\n")

		for comp in range(Minimum_decomposition,Maximum_decomposition):

			X = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(comp)+".npy")
			Y = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(comp+1)+".npy")

			#print("X", X.shape[1], "Y", Y.shape[1])

			for i in range(X.shape[1]):

				for j in range(Y.shape[1]):

					test = scipy.stats.pearsonr(X[:,i],Y[:,j])
					#print(i, j, test[0])
					outfile.write(str(comp)+"C_"+str(i+1)+"\t"+str(comp+1)+"C_"+str(j+1)+"\t"+str(abs(test[0]))+"\t")

					outfile.write(str(abs(Biton_sign_dict["M2_"+str(comp+1)+"C_"+str(j+1)][0]))+"\t"+str(abs(Biton_sign_dict["M3_"+str(comp+1)+"C_"+str(j+1)][0]))+"\t"+str(abs(Biton_sign_dict["M4_"+str(comp+1)+"C_"+str(j+1)][0]))+"\t"+str(abs(Biton_sign_dict["M5_"+str(comp+1)+"C_"+str(j+1)][0]))+"\t"+str(abs(Biton_sign_dict["M6_"+str(comp+1)+"C_"+str(j+1)][0]))+"\t"+str(abs(Biton_sign_dict["M7_"+str(comp+1)+"C_"+str(j+1)][0]))+"\t"+str(abs(Biton_sign_dict["M8_"+str(comp+1)+"C_"+str(j+1)][0]))+"\t"+str(abs(Biton_sign_dict["M9_"+str(comp+1)+"C_"+str(j+1)][0]))+"\t"+str(abs(Biton_sign_dict["M12_"+str(comp+1)+"C_"+str(j+1)][0]))+"\t"+str(abs(Biton_sign_dict["M13_"+str(comp+1)+"C_"+str(j+1)][0]))+"\t"+str(abs(Biton_sign_dict["M14_"+str(comp+1)+"C_"+str(j+1)][0]))+"\n")

	print("Job done!")

	return


def Load_complete_graph(datafolder, biton_folder, job, Minimum_decomposition, Maximum_decomposition):

	# Calculate the first nodes correlation with Biton
	Biton_M2_GC_CONTENT = pd.read_csv(biton_folder+'/Biton_M2_GC_CONTENT.txt', sep='\t')
	Biton_M3_SMOOTH_MUSCLE = pd.read_csv(biton_folder+'/Biton_M3_SMOOTH_MUSCLE.txt', sep='\t')
	Biton_M4_MITOCHONRIA_TRANSLATION = pd.read_csv(biton_folder+'/Biton_M4_MITOCHONRIA_TRANSLATION.txt', sep='\t')
	Biton_M5_INTERFERON = pd.read_csv(biton_folder+'/Biton_M5_INTERFERON.txt', sep='\t')
	Biton_M6_BASALLIKE = pd.read_csv(biton_folder+'/Biton_M6_BASALLIKE.txt', sep='\t')
	Biton_M7_CELLCYCLE = pd.read_csv(biton_folder+'/Biton_M7_CELLCYCLE.txt', sep='\t')
	Biton_M8_IMMUNE = pd.read_csv(biton_folder+'/Biton_M8_IMMUNE.txt', sep='\t')
	Biton_M9_UROTHELIALDIFF = pd.read_csv(biton_folder+'/Biton_M9_UROTHELIALDIFF.txt', sep='\t')
	Biton_M12_MYOFIBROBLASTS = pd.read_csv(biton_folder+'/Biton_M12_MYOFIBROBLASTS.txt', sep='\t')
	Biton_M13_BLCAPATHWAYS = pd.read_csv(biton_folder+'/Biton_M13_BLCAPATHWAYS.txt', sep='\t')
	Biton_M14_STRESS = pd.read_csv(biton_folder+'/Biton_M14_STRESS.txt', sep='\t')
	
	# Load data genes
	data = sc.read(datafolder+"/"+job+"_genes.txt",
							cache = False, first_column_names = True).transpose()
	#Biton_M2_GC_CONTENT
	sign_M2_idx = Biton_M2_GC_CONTENT.index.get_indexer(Biton_M2_GC_CONTENT.index.intersection(data.var.index,sort=None))
	plasma_M2_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M2_GC_CONTENT.index,sort=None))
	#Biton_M3_SMOOTH_MUSCLE
	sign_M3_idx = Biton_M3_SMOOTH_MUSCLE.index.get_indexer(Biton_M3_SMOOTH_MUSCLE.index.intersection(data.var.index,sort=None))
	plasma_M3_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M3_SMOOTH_MUSCLE.index,sort=None))
	#Biton_M4_MITOCHONRIA_TRANSLATION
	sign_M4_idx = Biton_M4_MITOCHONRIA_TRANSLATION.index.get_indexer(Biton_M4_MITOCHONRIA_TRANSLATION.index.intersection(data.var.index,sort=None))
	plasma_M4_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M4_MITOCHONRIA_TRANSLATION.index,sort=None))
	#Biton_M5_INTERFERON
	sign_M5_idx = Biton_M5_INTERFERON.index.get_indexer(Biton_M5_INTERFERON.index.intersection(data.var.index,sort=None))
	plasma_M5_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M5_INTERFERON.index,sort=None))
	#Biton_M6_BASALLIKE
	sign_M6_idx = Biton_M6_BASALLIKE.index.get_indexer(Biton_M6_BASALLIKE.index.intersection(data.var.index,sort=None))
	plasma_M6_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M6_BASALLIKE.index,sort=None))
	#Biton_M7_CELLCYCLE
	sign_M7_idx = Biton_M7_CELLCYCLE.index.get_indexer(Biton_M7_CELLCYCLE.index.intersection(data.var.index,sort=None))
	plasma_M7_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M7_CELLCYCLE.index,sort=None))
	#Biton_M8_IMMUNE
	sign_M8_idx = Biton_M8_IMMUNE.index.get_indexer(Biton_M8_IMMUNE.index.intersection(data.var.index,sort=None))
	plasma_M8_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M8_IMMUNE.index,sort=None))
	#Biton_M9_UROTHELIALDIFF
	sign_M9_idx = Biton_M9_UROTHELIALDIFF.index.get_indexer(Biton_M9_UROTHELIALDIFF.index.intersection(data.var.index,sort=None))
	plasma_M9_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M9_UROTHELIALDIFF.index,sort=None))
	#Biton_M12_MYOFIBROBLASTS
	sign_M12_idx = Biton_M12_MYOFIBROBLASTS.index.get_indexer(Biton_M12_MYOFIBROBLASTS.index.intersection(data.var.index,sort=None))
	plasma_M12_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M12_MYOFIBROBLASTS.index,sort=None))
	#Biton_M13_BLCAPATHWAYS
	sign_M13_idx = Biton_M13_BLCAPATHWAYS.index.get_indexer(Biton_M13_BLCAPATHWAYS.index.intersection(data.var.index,sort=None))
	plasma_M13_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M13_BLCAPATHWAYS.index,sort=None))
	#Biton_M14_STRESS
	sign_M14_idx = Biton_M14_STRESS.index.get_indexer(Biton_M14_STRESS.index.intersection(data.var.index,sort=None))
	plasma_M14_idx = data.var.index.get_indexer(data.var.index.intersection(Biton_M14_STRESS.index,sort=None))

	Biton_sign_dict = {}
	X = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(Minimum_decomposition)+".npy")
	
	for i in range(X.shape[1]):
		
		Biton_sign_dict["M2_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M2_idx],Biton_M2_GC_CONTENT["M2_GC_CONTENT"][sign_M2_idx])
		Biton_sign_dict["M3_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M3_idx],Biton_M3_SMOOTH_MUSCLE["M3_SMOOTH_MUSCLE"][sign_M3_idx])
		Biton_sign_dict["M4_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M4_idx],Biton_M4_MITOCHONRIA_TRANSLATION["M4_MITOCHONRIA_TRANSLATION"][sign_M4_idx])
		Biton_sign_dict["M5_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M5_idx],Biton_M5_INTERFERON["M5_INTERFERON"][sign_M5_idx])
		Biton_sign_dict["M6_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M6_idx],Biton_M6_BASALLIKE["M6_BASALLIKE"][sign_M6_idx])
		Biton_sign_dict["M7_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M7_idx],Biton_M7_CELLCYCLE["M7_CELLCYCLE"][sign_M7_idx])
		Biton_sign_dict["M8_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M8_idx],Biton_M8_IMMUNE["M8_IMMUNE"][sign_M8_idx])
		Biton_sign_dict["M9_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M9_idx],Biton_M9_UROTHELIALDIFF["M9_UROTHELIALDIFF"][sign_M9_idx])
		Biton_sign_dict["M12_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M12_idx],Biton_M12_MYOFIBROBLASTS["M12_MYOFIBROBLASTS"][sign_M12_idx])
		Biton_sign_dict["M13_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M13_idx],Biton_M13_BLCAPATHWAYS["M13_BLCAPATHWAYS"][sign_M13_idx])
		Biton_sign_dict["M14_"+str(Minimum_decomposition)+"C_"+str(i+1)] = scipy.stats.pearsonr(X[:,i][plasma_M14_idx],Biton_M14_STRESS["M14_STRESS"][sign_M14_idx])


	# Load the complete graph table
	Comp_corr_dict = {}
	Comp_sign_dict = {}
	with open('Graphs/'+job+"_full_ICA_graph.txt", 'r') as infile:
		first_line = infile.readline()
		for line in infile:
			tmp = line.split()
			# Get component-component correlations
			Comp_corr_dict[tmp[0]] = Comp_corr_dict.get(tmp[0],{})
			Comp_corr_dict[tmp[0]][tmp[1]] = float(tmp[2])
			# Get component-signature correlation
			Comp_sign_dict[tmp[1]] = Comp_sign_dict.get(tmp[1],tmp[3:])

	for i in range(Minimum_decomposition):
		Comp_sign_dict[str(Minimum_decomposition)+"C_"+str(i+1)] = [Biton_sign_dict["M2_"+str(Minimum_decomposition)+"C_"+str(i+1)][0],
																Biton_sign_dict["M3_"+str(Minimum_decomposition)+"C_"+str(i+1)][0],
																Biton_sign_dict["M4_"+str(Minimum_decomposition)+"C_"+str(i+1)][0],
																Biton_sign_dict["M5_"+str(Minimum_decomposition)+"C_"+str(i+1)][0],
																Biton_sign_dict["M6_"+str(Minimum_decomposition)+"C_"+str(i+1)][0],
																Biton_sign_dict["M7_"+str(Minimum_decomposition)+"C_"+str(i+1)][0],
																Biton_sign_dict["M8_"+str(Minimum_decomposition)+"C_"+str(i+1)][0],
																Biton_sign_dict["M9_"+str(Minimum_decomposition)+"C_"+str(i+1)][0],
																Biton_sign_dict["M12_"+str(Minimum_decomposition)+"C_"+str(i+1)][0],
																Biton_sign_dict["M13_"+str(Minimum_decomposition)+"C_"+str(i+1)][0],
																Biton_sign_dict["M14_"+str(Minimum_decomposition)+"C_"+str(i+1)][0]]

	# Get reverse component-component correlations
	Comp_corr_rev_dict = {}
	for key in Comp_corr_dict.keys():
		for key2 in Comp_corr_dict[key].keys():
			Comp_corr_rev_dict[key2] = Comp_corr_rev_dict.get(key2,{})
			Comp_corr_rev_dict[key2][key] = Comp_corr_dict[key][key2]

	return Comp_corr_dict, Comp_sign_dict, Comp_corr_rev_dict

def Get_K_ranked_descending_index(values_dict_value, k):
	
	return np.argsort(list(values_dict_value))[-k:][::-1]

def Test_gap(index_list, val_list, gap):
	
	kept_index_list = [index_list[0]]
	for i in range(len(index_list)-1):
		if (val_list[index_list[i]]/float(val_list[index_list[i+1]])) < gap:
			kept_index_list.append(index_list[i+1])
	
	return kept_index_list

def Get_index_to_comp_names(index_list, ref):
	
	tmp = []
	for val in index_list:
		tmp.append(str(ref)+"C_"+str(val+1))
	
	return tmp

def Perform_MNN(job, Minimum_decomposition, Maximum_decomposition, Gap, K, cor_thresh, Comp_corr_dict, Comp_sign_dict, Comp_corr_rev_dict):

	print("Applying MMN...")
	# Performing MNN with K neighbors
	MNN_dict = {}
	for order in range(Minimum_decomposition,Maximum_decomposition):
		for comp in range(order):
			IC_s = str(order)+"C_"+str(comp+1)

			#print("Testing neighbors of IC",IC_s)
			# Get K maximum value's index ranked descending or source-target
			rnk_s_list = Get_K_ranked_descending_index(Comp_corr_dict[IC_s].values(), K)
			#rnk_s_list = Test_gap(rnk_s_list, list(Comp_corr_dict[IC_s].values()), Gap)
			# Get K maximum value's index ranked descending of target-source
			IC_t_list = Get_index_to_comp_names(rnk_s_list, order+1)
			for IC_t in IC_t_list:
				rnk_t_list = Get_K_ranked_descending_index(Comp_corr_rev_dict[IC_t].values(), K)
				rnk_t_list = Test_gap(rnk_t_list, list(Comp_corr_rev_dict[IC_t].values()), Gap)
				IC_s_list = Get_index_to_comp_names(rnk_t_list, order)
				if IC_s in IC_s_list:
					MNN_dict[IC_s] = MNN_dict.get(IC_s,{})
					MNN_dict[IC_s][IC_t] = Comp_corr_dict[IC_s][IC_t]

	# Save MNN graph in file
	print("Saving full MMN graph to file...")
	with open('Graphs/'+job+"_MNN_full_graph.txt", 'w') as outfile:

		outfile.write("Source\tTarget\tCorrelation\tBiton_M2_GC_CONTENT\tBiton_M3_SMOOTH_MUSCLE\tBiton_M4_MITOCHONRIA_TRANSLATION\tBiton_M5_INTERFERON\tBiton_M6_BASALLIKE\tBiton_M7_CELLCYCLE\tBiton_M8_IMMUNE\tBiton_M9_UROTHELIALDIFF\tBiton_M12_MYOFIBROBLASTS\tBiton_M13_BLCAPATHWAYS\tBiton_M14_STRESS\n")

		for source in MNN_dict.keys():
			for target in MNN_dict[source].keys():
				if MNN_dict[source][target] > cor_thresh:
					outfile.write(source+"\t"+target+"\t"+str(MNN_dict[source][target]))
					for corr in Comp_sign_dict[target]:
						outfile.write("\t"+str(corr))
					outfile.write("\n")

	# Load the MNN graph table
	MNN_corr_dict = {}
	MNN_sign_dict = {}
	with open('Graphs/'+job+"_MNN_full_graph.txt", 'r') as infile:
		first_line = infile.readline()
		for line in infile:
			tmp = line.split()
			# Get component-component correlations
			MNN_corr_dict[tmp[0]] = MNN_corr_dict.get(tmp[0],{})
			MNN_corr_dict[tmp[0]][tmp[1]] = float(tmp[2])
			# Get component-signature correlation
			MNN_sign_dict[tmp[1]] = MNN_sign_dict.get(tmp[1],tmp[3:])

	return MNN_corr_dict, MNN_sign_dict

def Convert_MMN_dict_to_NetworkX(MNN_corr_dict):
	
	print("Converting dict to network...")
	### Generate NetworkX directed graph from MNN
	Directed_graph=nx.DiGraph()
	# Adding Nodes
	for source in MNN_corr_dict.keys():
		Directed_graph.add_node(source)
		for target in MNN_corr_dict[source].keys():
			Directed_graph.add_node(target)
			Directed_graph.add_edge(source,target,weight=MNN_corr_dict[source][target])

	return Directed_graph

def Delete_splits(Directed_graph):

	#Copy graph for reference and to use as output without changing the input one
	DG_nosplit = Directed_graph.copy()
	### Find nodes with out/in degree superior to 1 and delete outgoing/ingoing edges (aka splits)

	# Delete outgoing degree edges
	print("Deleting outgoing edges...")
	remove_nodes = [node for node,degree in dict(Directed_graph.out_degree()).items() if degree > 1]
	#print(remove_nodes)
	remove_edges = Directed_graph.out_edges(remove_nodes)
	#print(remove_edges)
	DG_nosplit.remove_edges_from(remove_edges)

	# Delete ingoing degree edges
	print("Deleting ingoing edges...")
	remove_nodes = [node for node,degree in dict(Directed_graph.in_degree()).items() if degree > 1]
	#print(remove_nodes)
	remove_edges = Directed_graph.in_edges(remove_nodes)
	#print(remove_edges)
	DG_nosplit.remove_edges_from(remove_edges)

	return DG_nosplit

def Clean_non_persistent_components(DG_nosplit, len_thresh=10):
	
	print("Prunning non persistent components...")
	# Check for connected components and remove based on a threshold
	#CC = nx.connected_components(G_nosplit)
	CC = nx.weakly_connected_components(DG_nosplit)
	remove_cc_nodes = []
	count = 0
	i = 0
	for c in CC:
		if len(c) < len_thresh:
			remove_cc_nodes.extend(list(c))
			count += 1
		i += 1

	# Store the cleaned graph
	DG_persist = DG_nosplit.copy()
	DG_persist.remove_nodes_from(remove_cc_nodes)

	return DG_persist

### Put nodes attributes in dictionary of dictionary
def Add_node_attributes(graph, Comp_sign_dict):
	
	Nodes_dict = {}
	Sign_annot = ("Biton_M2_GC_CONTENT","Biton_M3_SMOOTH_MUSCLE","Biton_M4_MITOCHONRIA_TRANSLATION","Biton_M5_INTERFERON","Biton_M6_BASALLIKE","Biton_M7_CELLCYCLE","Biton_M8_IMMUNE","Biton_M9_UROTHELIALDIFF","Biton_M12_MYOFIBROBLASTS","Biton_M13_BLCAPATHWAYS","Biton_M14_STRESS")
	for node in list(graph.nodes()):
		Nodes_dict[node] = {}
		for annot in range(len(Comp_sign_dict[node])):
			Nodes_dict[node][Sign_annot[annot]] = float(Comp_sign_dict[node][annot])
	nx.set_node_attributes(graph, Nodes_dict)
	
	return graph

### Put Edges attributes in dictionary of dictionary
def Add_edge_attributes(graph, job):
	
	Edges_dict = {}
	Sign_annot = ("Biton_M2_GC_CONTENT","Biton_M3_SMOOTH_MUSCLE","Biton_M4_MITOCHONRIA_TRANSLATION","Biton_M5_INTERFERON","Biton_M6_BASALLIKE","Biton_M7_CELLCYCLE","Biton_M8_IMMUNE","Biton_M9_UROTHELIALDIFF","Biton_M12_MYOFIBROBLASTS","Biton_M13_BLCAPATHWAYS","Biton_M14_STRESS")
	for edge in list(graph.edges()):
		Edges_dict[edge] = {}
		start_order = int(edge[0].split("C_")[0])
		start_ic_num = int(edge[0].split("C_")[-1])
		Cs = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(start_order)+".npy")
		end_order = int(edge[1].split("C_")[0])
		end_ic_num = int(edge[1].split("C_")[-1])
		Ce = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(end_order)+".npy")
		s_e_corr = scipy.stats.pearsonr(Cs[:,start_ic_num-1],Ce[:,end_ic_num-1])
		Edges_dict[edge]["Pear_cor"] = abs(s_e_corr[0])
		Edges_dict[edge]["Length"] = end_order-start_order
		Edges_dict[edge]["Average_metagene"] = (Cs[:,start_ic_num-1] + Ce[:,end_ic_num-1])/2

		for attrib in Sign_annot:
			Edges_dict[edge]["Average_"+attrib] = (graph.nodes[edge[0]][attrib]+graph.nodes[edge[1]][attrib])/2
	nx.set_edge_attributes(graph, Edges_dict)
	
	return graph


def Simplify_persistent_components(DG_persist, Comp_sign_dict, job):

	print("Simplifying persistent components...")
	# Find begining and ending of each graph connected components
	CC_dict = {} # Store straight lines info: start IC, end IC, line length
	CC_dict_start = {} # Dict with key=start node and value=end node
	CC_dict_end = {}# Dict with key=end node and value=start node
	Start_nodes = []
	End_nodes = []
	CC = nx.weakly_connected_components(DG_persist)
	for ccomp in CC:
		for node in ccomp:
			if DG_persist.in_degree(node) < 1:
				start = node
				Start_nodes.append(node)
			elif DG_persist.out_degree(node) < 1:
				end = node
				End_nodes.append(node)
		CC_dict[(start,end)] = len(ccomp)
		CC_dict_start[start] = end
		CC_dict_end[end] = start


	### Simplify components by their start and end nodes
	DG_persist_simple = nx.Graph()
	# Sort connected component by lowest order first
	Start_nodes.sort(key = lambda x: int(x.split("C_")[0]))
	for s_node in Start_nodes:
		e_node = CC_dict_start[s_node]

		DG_persist_simple.add_node(s_node)
		DG_persist_simple.add_node(e_node)
		DG_persist_simple.add_edge(s_node,e_node)

	
	DG_persist_simple = Add_node_attributes(DG_persist_simple, Comp_sign_dict)
	DG_persist_simple = Add_edge_attributes(DG_persist_simple, job)

	return DG_persist_simple, CC_dict_start

def Recontruct_splits(DG_persist_simple, job, Comp_sign_dict, CC_dict_start, Corr_threshold=0.3, Gap=1.5, K=0, Penalty=0.05):

	print("Reconstructing splits...")
	DG_persist_rewired = DG_persist_simple.copy()

	# Calculate correlations between start/end nodes and create ALL edges above threshold and Gap
	test_count = 0
	overlap_count = 0
	no_candidate_count = 0
	gap_count = 0
	add_count = 0
	Added_edges = []
	Added_edges_cor = []

	start_nodes_list = list(CC_dict_start.keys())
	for start_1 in start_nodes_list:
		C1_end = CC_dict_start[start_1]
		end_order = int(C1_end.split("C_")[0])
		end_ic_num = int(C1_end.split("C_")[-1])
		C1 = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(end_order)+".npy")
		candidate_link_list = []
		candidate_cor_list = []
		for start_2 in start_nodes_list:
			# Avoid creating loops in same component
			if start_1 != start_2:
				C2_start = start_2
				start_order = int(C2_start.split("C_")[0])
				start_ic_num = int(C2_start.split("C_")[-1])
				# Avoid loops in overlapping components
				if end_order < start_order:
					C2 = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(start_order)+".npy")
					# Test correlation
					test_cor = scipy.stats.pearsonr(C1[:,end_ic_num-1],C2[:,start_ic_num-1])
					# Correct correlation to prioritise close links by putting a penalty
					alpha = Penalty
					penalty_cor = abs(test_cor[0]) - alpha * (start_order-end_order)
					#if abs(test_cor[0]) > Corr_threshold: # Before penalty
					if penalty_cor > Corr_threshold:
						#print(C2_start, C1_end, abs(test_cor[0]), penalty_cor)
						candidate_link_list.append(C2_start)
						#candidate_cor_list.append(abs(test_cor[0])) #Before penalty
						candidate_cor_list.append(abs(penalty_cor))
						#G_cleaned_rewired.add_edge(C1_end,C2_start,weight=abs(test_cor[0]))
				else:
					overlap_count += 1
			else:
				test_count+=1
		# Check if multiple candidate links are found
		if len(candidate_link_list)>=2:
			# Take only K maximum correlation edges
			#print(candidate_link_list)
			#print(candidate_cor_list)
			### rnk_link_ind = Get_K_ranked_descending_index(candidate_cor_list,K)
			if K:
				rnk_link_ind = Get_K_ranked_descending_index(candidate_cor_list,K)
			else:
				rnk_link_ind = Get_K_ranked_descending_index(candidate_cor_list,len(candidate_cor_list))

			# Take the maximum correlation edge always at first
			DG_persist_rewired.add_edge(C1_end,candidate_link_list[0],weight=candidate_cor_list[0])
			add_count += 1
			Added_edges.append((C1_end,candidate_link_list[0]))
			Added_edges_cor.append(candidate_cor_list[0])
			for cand_ind in range(len(rnk_link_ind)-1):
				# Test if gap between other candidate edges is not beyond threshold
				if candidate_cor_list[rnk_link_ind[0]]/candidate_cor_list[rnk_link_ind[cand_ind+1]] < Gap:
					DG_persist_rewired.add_edge(C1_end,candidate_link_list[rnk_link_ind[cand_ind+1]],weight=candidate_cor_list[rnk_link_ind[cand_ind+1]])
					add_count += 1
					Added_edges.append((C1_end,candidate_link_list[rnk_link_ind[cand_ind+1]]))
					Added_edges_cor.append(candidate_cor_list[rnk_link_ind[cand_ind+1]])
				else:
					gap_count += 1
		# Simply keep the only candidate link found
		elif candidate_link_list:
			DG_persist_rewired.add_edge(C1_end,candidate_link_list[0],weight=candidate_cor_list[0])
			add_count += 1
			Added_edges.append((C1_end,candidate_link_list[0]))
			Added_edges_cor.append(candidate_cor_list[0])
		else:
			no_candidate_count+=1

	DG_persist_rewired = Add_node_attributes(DG_persist_rewired, Comp_sign_dict)
	DG_persist_rewired = Add_edge_attributes(DG_persist_rewired, job)
	print("Saving processed graph in file...")
	if K:
		nx.write_weighted_edgelist(DG_persist_rewired, 'Graphs/'+job+'_MNN_cleaned_graph_'+str(K)+'_edges_only.txt')
	else:
		nx.write_weighted_edgelist(DG_persist_rewired, 'Graphs/'+job+'_MNN_cleaned_graph_ALL_edges.txt')
	
	#print("Same component:", test_count)
	#print("Overlaps:", overlap_count)
	#print("No Candidates:", no_candidate_count)
	#print("Gaps:", gap_count)
	#print("Additions:", add_count)
	print("Job done!")
	
	return DG_persist_rewired, Added_edges
