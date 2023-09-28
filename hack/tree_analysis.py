import pandas as pd
import numpy as np
import networkx as nx
import scipy.stats
import matplotlib.pyplot as plt
from sklearn import preprocessing

def Alternating_left_right_node_coord(graph, comp_start_dict):
	
	start_nodes_list = list(comp_start_dict.keys())
	start_nodes_list.sort(key = lambda x: int(x.split("C_")[0]))
	Nodes_coord_dict = {}
	x_coord = 1	
	mult_switch = 1
	for s_node in start_nodes_list:
		e_node = comp_start_dict[s_node]
		x = (x_coord//2)*mult_switch
		x_coord += 1
		mult_switch *= -1

		Nodes_coord_dict[s_node] = {}
		Nodes_coord_dict[s_node]["pos"] = (x,int(s_node.split("C_")[0]))
		Nodes_coord_dict[e_node] = {}
		Nodes_coord_dict[e_node]["pos"] = (x,int(e_node.split("C_")[0]))
	nx.set_node_attributes(graph, Nodes_coord_dict)
	
	return graph

def Slide_right_node_coord(graph, comp_start_dict):
	
	start_nodes_list = list(comp_start_dict.keys())
	start_nodes_list.sort(key = lambda x: int(x.split("C_")[0]))
	Nodes_coord_dict = {}
	x_coord = 1
	for s_node in start_nodes_list:
		e_node = comp_start_dict[s_node]
		x = x_coord
		x_coord += 1

		Nodes_coord_dict[s_node] = {}
		Nodes_coord_dict[s_node]["pos"] = (x,int(s_node.split("C_")[0]))
		Nodes_coord_dict[e_node] = {}
		Nodes_coord_dict[e_node]["pos"] = (x,int(e_node.split("C_")[0]))
	nx.set_node_attributes(graph, Nodes_coord_dict)
	
	return graph

def Get_signature_correlations(graph, component_str, signature_choice, gap_threshold):
	
	# Plot signature correlations barplot
	signatures=["Biton_M2_GC_CONTENT","Biton_M3_SMOOTH_MUSCLE","Biton_M4_MITOCHONRIA_TRANSLATION","Biton_M5_INTERFERON","Biton_M6_BASALLIKE","Biton_M7_CELLCYCLE","Biton_M8_IMMUNE","Biton_M9_UROTHELIALDIFF","Biton_M12_MYOFIBROBLASTS","Biton_M13_BLCAPATHWAYS","Biton_M14_STRESS"]
	sign_ind = signatures.index(signature_choice)
	score_list = []
	for sign in signatures:
		score_list.append(graph.nodes[component_str][sign])
		#print(graph.nodes[component_str][sign])
	sign_score = score_list[sign_ind]
	score_list.sort(reverse=True)
	print(score_list, sign_score)
	if sign_score == score_list[0]:
		if sign_score/score_list[1] > gap_threshold:
			return True
		else:
			return False
	else:
		if score_list[0]/sign_score > gap_threshold:
			return True
		else:
			return False
	
def Get_SD_genes(component_str, SD_threshold, foldername, job):
	
	gene_list = []
	Genes = pd.read_csv(foldername+"/"+job+"_genes.txt", sep='\t')
	order = int(component_str.split("C_")[0])
	IC = int(component_str.split("C_")[1])
	S = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(order)+".npy")
	
	IC_vect = S[:,IC-1]
	# Extract genes
	IC_std = np.std(IC_vect)
	IC_mean = np.mean(IC_vect)
	#print(IC_std)
	# Genes above SD
	IC_plus = IC_vect > IC_mean+(IC_std*SD_threshold)
	#print(Genes[IC_plus].values[0])
	gene_list.extend(Genes[IC_plus]["Genes"].to_list())
	# Genes below SD
	IC_minus = IC_vect < IC_mean-(IC_std*SD_threshold)
	gene_list.extend(Genes[IC_minus]["Genes"].to_list())
	
	return gene_list

def Create_plot_colors(graph, max_dim, layout, threshold_sign, signature, edge_select, print_node_labels,
						print_edge_labels, print_all_nodes, source_node, node_cmap, edge_cmap,
						genes_test, SD, genes_prop, file_name, foldername, job, added_edges, start_nodes,show):
	
	# Find nodes and edges if 1 node has high correlation with signature
	edge_width = 3
	threshold_sign = float(threshold_sign)
	SD = int(SD)
	genes_prop = float(genes_prop)
	nodes_to_draw = []
	edges_to_draw = []
	nodes_gene_prop = []
	# Draw nodes based on signature correlation
	if len(source_node)<1:
		# Find sure nodes
		for node in graph.nodes():
			if graph.nodes()[node][signature] > threshold_sign:
				# To look later to select sign nodes based on relative gap of sign corr
				"""
				#Gap_sign_test = Get_signature_correlations(graph, node, signature, 1.5)
				#print(Gap_sign_test, node)
				"""
				nodes_to_draw.append(node)
		# Get all connected nodes
		if print_all_nodes:
			for node in nodes_to_draw:
				sign_nodes = nx.node_connected_component(graph, node)
				for s_node in sign_nodes:
					if s_node not in nodes_to_draw:
						nodes_to_draw.append(s_node)
		# Get edges containing nodes to draw
		for edge in graph.edges():
			if edge[0] in nodes_to_draw:
				edges_to_draw.append(edge)
			elif edge[1] in nodes_to_draw:
				edges_to_draw.append(edge)
	# Draw nodes from input list
	else:
		nodes_to_draw = source_node.split(",")
		# Get all connected nodes
		for node in nodes_to_draw:
			sign_nodes = nx.node_connected_component(graph, node)
			for s_node in sign_nodes:
				if s_node not in nodes_to_draw:
					nodes_to_draw.append(s_node)
		# Get edges containing nodes to draw
		for edge in graph.edges():
			if edge[0] in nodes_to_draw:
				edges_to_draw.append(edge)
			elif edge[1] in nodes_to_draw:
				edges_to_draw.append(edge)
	# Get nodes containing genes of interest
	if len(genes_test)>0:
		gene_list = genes_test.split(",")
		print("Looking for genes of interest")
		nodes_to_draw = []
		for node in graph.nodes():
			SD_gene_list = Get_SD_genes(node, SD, foldername, job)
			gene_count = 0
			for gene in gene_list:
				if gene in SD_gene_list:
					gene_count += 1
			if gene_count>0:
				prop_found = gene_count/len(gene_list)
				if prop_found >= genes_prop:
					nodes_to_draw.append(node)
					nodes_gene_prop.append(prop_found)
		# Get edges containing nodes to draw
		edges_to_draw = []
		for edge in graph.edges():
			if edge[0] in nodes_to_draw:
				edges_to_draw.append(edge)
			elif edge[1] in nodes_to_draw:
				edges_to_draw.append(edge)
		#print("Found nodes:", len(nodes_to_draw))
		#print("Found edges:", len(edges_to_draw))
	
	if len(nodes_to_draw) < 1:
		print("No components found...")
		return
	else:
		print("Found", len(nodes_to_draw), "Components!")
	
	# Get nodes labels
	node_labels = {}
	for node in nodes_to_draw:
		node_labels[node] = node
	# Get nodes color based on signature correlation
	node_colors = []
	for node in nodes_to_draw:
		node_colors.append(graph.nodes()[node][signature])
	# Get edges width based on average signature correlation
	edge_widths = []
	for edge in edges_to_draw:
		edge_widths.append(graph.edges()[edge]["Average_"+signature])
	# Get edges color based on average signature correlation
	edge_colors = []
	edge_labels_sign = {}
	for edge in edges_to_draw:
		edge_colors.append(graph.edges()[edge]["Average_"+signature])
		edge_labels_sign[edge] = round(graph.edges()[edge]["Average_"+signature], 2)
	# Get edges color based on Pearson correlation
	edge_colors_cor = []
	edge_labels_pear = {}
	for edge in edges_to_draw:
		edge_colors_cor.append(graph.edges()[edge]["Pear_cor"])
		edge_labels_pear[edge] = round(graph.edges()[edge]["Pear_cor"], 2)
	# Get edges style: solid if true component and dashed if added
	edge_dashed_name = []
	edge_solid_name = []
	edge_dashed_colors = []
	edge_solid_colors = []
	edge_dashed_colors_cor = []
	edge_solid_colors_cor = []
	edge_dashed_labels_sign = {}
	edge_solid_labels_sign = {}
	edge_dashed_labels_pear = {}
	edge_solid_labels_pear = {}
	count = 0
	for edge in edges_to_draw:
		if edge in added_edges:
			edge_dashed_name.append(edge)
			edge_dashed_colors.append(edge_colors[count])
			edge_dashed_colors_cor.append(edge_colors_cor[count])
			edge_dashed_labels_sign[edge] = edge_labels_sign[edge]
			edge_dashed_labels_pear[edge] = edge_labels_pear[edge]
		else:
			edge_solid_name.append(edge)
			edge_solid_colors.append(edge_colors[count])
			edge_solid_colors_cor.append(edge_colors_cor[count])
			edge_solid_labels_sign[edge] = edge_labels_sign[edge]
			edge_solid_labels_pear[edge] = edge_labels_pear[edge]
		count+=1
	
	# Plot the selected components
	if layout == "Alternating":
		graph = Alternating_left_right_node_coord(graph, start_nodes)
	if layout == "Climbing":
		graph = Slide_right_node_coord(graph, start_nodes)
	
	# Get positions of nodes and consequently edges
	pos=nx.get_node_attributes(graph,'pos')
	# Get node labels positions on top-right
	pos_labels={}
	for p in pos:
		pos_labels[p] = (pos[p][0]+3,pos[p][1]+1)
	
	fig, ax = plt.subplots()
	
	if print_node_labels:
		nx.draw_networkx_labels(graph,pos,labels=node_labels,
											font_size=10,font_weight='bold',ax=ax)
		nodes_fig = nx.draw_networkx_nodes(graph, pos,
						nodelist=nodes_to_draw, node_size=0,
						node_color=node_colors, cmap=plt.cm.get_cmap(node_cmap), vmin=0, vmax=1,
						edgecolors="k",linewidths=2,ax=ax)
		alpha_val = 0.5
	elif len(genes_test)>0:
		nodes_fig = nx.draw_networkx_nodes(graph, pos,
						nodelist=nodes_to_draw, node_size=100,
						node_color=nodes_gene_prop, cmap=plt.cm.get_cmap(node_cmap), vmin=0, vmax=1,
						edgecolors="k",linewidths=2,ax=ax)
		alpha_val = 1
	else:
		nodes_fig = nx.draw_networkx_nodes(graph, pos,
						nodelist=nodes_to_draw, node_size=100,
						node_color=node_colors, cmap=plt.cm.get_cmap(node_cmap), vmin=0, vmax=1,
						edgecolors="k",linewidths=2,ax=ax)
		alpha_val = 1

	### Test for dashed edges
	if edge_select == "Signature_correlation":
		# Solid edges first
		if edge_solid_name:
			edges_border = nx.draw_networkx_edges(graph, pos, edgelist=edge_solid_name,
								  width=edge_width*1.5, edge_color="k", style='solid',alpha=alpha_val)
			edges_fig = nx.draw_networkx_edges(graph, pos, edgelist=edge_solid_name, width=edge_width,
								   edge_color=edge_solid_colors, style='solid',alpha=alpha_val,
								   edge_cmap=plt.cm.get_cmap(edge_cmap), edge_vmin=0, edge_vmax=1)
		# Dashed edges second
		if edge_dashed_name:
			edges_border = nx.draw_networkx_edges(graph, pos, edgelist=edge_dashed_name,
								  width=edge_width*0.5, edge_color="k", style='solid',alpha=alpha_val)
			edges_fig = nx.draw_networkx_edges(graph, pos, edgelist=edge_dashed_name, width=edge_width*0.5,
								   edge_color=edge_dashed_colors, style='dashed',alpha=alpha_val,
								   edge_cmap=plt.cm.get_cmap(edge_cmap), edge_vmin=0, edge_vmax=1)
		if print_edge_labels:
			nx.draw_networkx_edge_labels(graph, pos,
										edge_labels=edge_labels_sign,label_pos=0.5,
										font_weight='bold')
		
	elif edge_select == "Pearson_correlation":
		# Solid edges first
		if edge_solid_name:
			edges_border = nx.draw_networkx_edges(graph, pos, edgelist=edge_solid_name,
								  width=edge_width*1.5, edge_color="k", style='solid',alpha=alpha_val)
			edges_fig = nx.draw_networkx_edges(graph, pos, edgelist=edge_solid_name, width=edge_width,
								   edge_color=edge_solid_colors_cor, style='solid',alpha=alpha_val,
								   edge_cmap=plt.cm.get_cmap(edge_cmap), edge_vmin=0, edge_vmax=1)
		# Dashed edges second
		if edge_dashed_name:
			edges_border = nx.draw_networkx_edges(graph, pos, edgelist=edge_dashed_name,
								  width=edge_width*0.5, edge_color="w", style='solid',alpha=alpha_val)
			edges_fig = nx.draw_networkx_edges(graph, pos, edgelist=edge_dashed_name, width=edge_width*0.5,
								   edge_color=edge_dashed_colors_cor, style='dashed',alpha=alpha_val,
								   edge_cmap=plt.cm.get_cmap(edge_cmap), edge_vmin=0, edge_vmax=1)
		if print_edge_labels:
			nx.draw_networkx_edge_labels(graph, pos,
										edge_labels=edge_labels_pear,label_pos=0.5,
										font_weight='bold')
	### End of test for dashed edges
	
	
	plt.ylabel('Decomposition Order',weight='bold',fontsize=18)
	plt.yticks(np.arange(0, max_dim+5, 5),labels=np.arange(0, max_dim+5, 5), fontsize=14, fontweight='bold')
	plt.grid(axis='y', linestyle='--')
	plt.grid(visible=None,axis='x')
	#plt.sci(edges_border)    
	plt.sci(edges_fig)
	if edge_select == "Signature_correlation":
		plt.colorbar().set_label("Edges Average "+signature+" correlation",weight='bold',fontsize=18)
	elif edge_select == "Pearson_correlation":
		plt.colorbar().set_label("Edges Pearson correlation",weight='bold',fontsize=18)
	plt.gci().set_cmap(edge_cmap)
	plt.sci(nodes_fig)
	if len(genes_test)>0:
		plt.colorbar().set_label("Nodes found genes proportion",weight='bold',fontsize=18)
	else:
		plt.colorbar().set_label("Nodes "+signature+" correlation",weight='bold',fontsize=18)
	plt.gci().set_cmap(node_cmap)
	ax.tick_params(left=True, bottom=False, labelleft=True, labelbottom=False)
	#ax.collections[1].set_edgecolor("#FF0000")
	#ax.collections[0].set_edgecolor("#FF0000") 
	plt.savefig(file_name+".svg")
	if show:
		plt.show()
	
	return

def Orient_component(comp_vect, score_type, threshold, orient=False):
	
	if score_type == "zscore":
		comp_vect_zscr = scipy.stats.zscore(comp_vect)
		if orient:
			top_thresh = sum([x for x in comp_vect_zscr if x > threshold])
			bot_thresh = sum([x for x in comp_vect_zscr if x < -threshold])
			if abs(top_thresh) > abs(bot_thresh):
				return comp_vect_zscr
			else:
				return comp_vect_zscr * -1
		else:
			return comp_vect_zscr
	elif score_type == "icascore":
		comp_vect_sort = sorted(comp_vect)
		if orient:
			top_thresh = sum(comp_vect_sort[-threshold:])
			bot_thresh = sum(comp_vect_sort[:threshold])
			if abs(top_thresh) > abs(bot_thresh):
				return comp_vect
			else:
				return comp_vect * -1
		else:
			return comp_vect
	elif score_type == "rank":
		comp_vect_zscr = scipy.stats.zscore(comp_vect)
		if orient:
			top_thresh = sum([x for x in comp_vect_zscr if x > threshold])
			bot_thresh = sum([x for x in comp_vect_zscr if x < -threshold])
			if abs(top_thresh) > abs(bot_thresh):
				comp_vect = comp_vect * -1
				tmp = comp_vect.argsort()
				ranks = np.empty_like(tmp)
				ranks[tmp] = np.arange(len(comp_vect))
				return ranks + 1
			else:
				tmp = comp_vect.argsort()
				ranks = np.empty_like(tmp)
				ranks[tmp] = np.arange(len(comp_vect))
				return ranks + 1
		else:
			tmp = comp_vect.argsort()
			ranks = np.empty_like(tmp)
			ranks[tmp] = np.arange(len(comp_vect))
			return ranks + 1
	else:
		print("Wrong score type input... Input untouched!")
		return comp_vect

def Compute_average_segments(graph_full, graph_simple, score_type, threshold, datafolder, job, orient=False):
	
	# Load gene names
	genes_pd = pd.read_csv(datafolder+"/"+job+"_genes.txt", sep='\t')
	# Extract full component paths from NetworkX graph
	graph_full_fcp = nx.weakly_connected_components(graph_full)
	graph_simple_fcp = nx.connected_components(graph_simple)
	# Loop over segments
	print("Computing average scores...")
	seg_average_list = []
	for segment in graph_full_fcp:
		seg_all_vect = []
		# Loop over components
		for component in segment:
			tmp = component.split("C_")
			comp = int(tmp[0])
			IC = int(tmp[1])
			X = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(comp)+".npy")
			IC_vect = X[:,IC-1]
			seg_all_vect.append(IC_vect)
		# Compute average
		sum_vect = np.zeros(len(seg_all_vect[0]))
		for vect in seg_all_vect:
			# Orient the component based on score type and threshold
			orient_average_vect = Orient_component(vect, score_type, threshold, orient)
			# Add the oriented score into the list
			sum_vect += orient_average_vect
		# Average the list of scores
		average_sum_vect = sum_vect/(len(segment))
		
		if score_type == "rank":
			average_sum_vect = np.rint(average_sum_vect)
		
		seg_average_list.append(average_sum_vect)
	# Get segment simplified names
	segment_simple_names_list = []
	for segment in graph_simple_fcp:
		comp_pair = []
		for comp in segment:
			comp_pair.append(comp)
		comp1 = comp_pair[0].split("C_")[0]
		comp2 = comp_pair[1].split("C_")[0]
		if int(comp1) < int(comp2):
			segment_simple_names_list.append(comp_pair[0]+"to"+comp_pair[1])
		else:
			segment_simple_names_list.append(comp_pair[1]+"to"+comp_pair[0])

	print("Saving results...")
	# Store the components in a pandas table for export
	average_segment_pd = pd.DataFrame.from_dict(dict(zip(segment_simple_names_list, seg_average_list)))
	#print(average_segment_pd[segment_simple_names_list[0]])
	average_segment_pd.index = genes_pd["Genes"].to_list()

	average_segment_pd.to_csv('Graphs/'+job+"_"+score_type+"_thresh"+str(threshold)+"_average_segments.txt", sep="\t", index_label="Genes")
	print("Job done!")
	
	return

def Compute_simplified_normalised_ranks(persist_graph, datafolder, job, threshold=0.0):

	genes_pd = pd.read_csv(datafolder+"/"+job+"_genes.txt", sep='\t')

	print("Generating stable ranks...")
	Detailed_persistent_comp = nx.weakly_connected_components(persist_graph)
	Simplified_comp_dict = {}
	for persistent_comp in Detailed_persistent_comp:
		# Get simplified persistent component name
		dim_list = []
		comp_list = []
		for comp in persistent_comp:
			comp_list.append(comp)
			dim_list.append(int(comp.split("C_")[0]))
		dim_min = min(dim_list)
		dim_max = max(dim_list)
		simp_comp = str(comp_list[dim_list.index(dim_min)])+"_to_"+str(comp_list[dim_list.index(dim_max)])
		#print(simp_comp)

		All_comp_array = []
		for comp in persistent_comp:
			tmp = comp.split("C_")
			order = tmp[0]
			N = tmp[1]
			C = np.load(job+"_deconvolutions/"+job+"_S_matrix_C"+str(order)+".npy")
			ica_score_comp = C[:,int(N)-1]
			All_comp_array.append(ica_score_comp)
			#print(ica_score_comp)
		All_comp_array = np.array(All_comp_array)

		# Get ranks of components in ascending order: most negative = 1, most positive = max_len
		All_comp_array_rank = []
		for comp in range(len(All_comp_array)):
			tmp = All_comp_array[comp].argsort()
			ranks = np.empty_like(tmp)
			ranks[tmp] = np.arange(len(All_comp_array[comp]))
			ranks = ranks + 1
			All_comp_array_rank.append(ranks)

		All_comp_array_rank = np.array(All_comp_array_rank)

		All_comp_array_rank_norm = preprocessing.minmax_scale(All_comp_array_rank, feature_range=(-1, 1), axis=1)
		# Simplify normalised ranks to -1, 0 or +1 based on the set threshold
		All_comp_array_rank_norm_clean = np.where(All_comp_array_rank_norm >= -threshold, All_comp_array_rank_norm, -1)
		All_comp_array_rank_norm_clean = np.where(All_comp_array_rank_norm_clean <= threshold, All_comp_array_rank_norm_clean, 1)
		All_comp_array_rank_norm_clean = np.where((All_comp_array_rank_norm_clean > -threshold) ^ (All_comp_array_rank_norm_clean < threshold), All_comp_array_rank_norm_clean, 0)
		# Find genes that are not only -1 or +1
		all_neg = np.all(All_comp_array_rank_norm_clean == -1, axis=0)
		all_pos = np.all(All_comp_array_rank_norm_clean == 1, axis=0)
		not_interest = ~all_neg & ~all_pos
		# Calculate mean normalised rank for only positive or negative and put 0 for rest
		Mean_norm_rank_cleaned = All_comp_array_rank_norm.mean(0)
		# Put all non stable values to zero
		Mean_norm_rank_cleaned[not_interest] = int(0)
		Simplified_comp_dict[simp_comp] = Mean_norm_rank_cleaned

	print("Saving file...")
	Simplified_norm_ranks_df = pd.DataFrame.from_dict(Simplified_comp_dict)
	Simplified_norm_ranks_df.index = genes_pd["Genes"].to_list()
	Simplified_norm_ranks_df.to_csv('Graphs/'+job+"_thresh_"+str(threshold)+"_average_stability_persistent_components.txt",
									sep="\t", index_label="Genes", float_format='%.5f')
	print("Done!")
	return
