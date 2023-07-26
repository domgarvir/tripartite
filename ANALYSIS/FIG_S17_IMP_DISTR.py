#!/usr/bin/env python3
from Mnetworks import *
from Metrics import *
from Graphics import *

#First get ranking of all plants with their actual value of correlation
#we want to plot the distribution of importances

filtered_network_names = ['Sinohara_1_ALL_PH', 'Sinohara_2_ALL_PH', 'Sinohara_3_ALL_PH', 'Sinohara_4_ALL_PH', 'Sinohara_ALL_A_PH', 'Sinohara_ALL_E_PH', 'Sinohara_ALL_I_PH', 'Sinohara_2_E_PH','Sinohara_3_E_PH', 'Sinohara_4_I_PH', 'Melian_OO_OO_PH', 'Hackett_1_ALL_PH', 'Hackett_2_ALL_PH', 'Hackett_1_S_PH', 'Hackett_1_GL_PH', 'Pocock_OO_OO_PH','Melian_OO_OO_HSD', 'McFayden_ALL_A_HPa', 'McFayden_1_A_HPa', 'McFayden_2_A_HPa','McFayden_3_A_HPa', 'McFayden_4_A_HPa', 'McFayden_5_A_HPa', 'McFayden_6_A_HPa', 'McFayden_7_A_HPa', 'McFayden_8_A_HPa', 'McFayden_9_A_HPa', 'McFayden_10_A_HPa', 'McFayden_ALL_B_HPa', 'McFayden_1_B_HPa', 'McFayden_2_B_HPa', 'McFayden_3_B_HPa',  'McFayden_4_B_HPa', 'McFayden_5_B_HPa', 'McFayden_6_B_HPa', 'McFayden_7_B_HPa','McFayden_8_B_HPa', 'McFayden_9_B_HPa', 'McFayden_10_B_HPa', 'Hackett_1_ALL_HPa','Hackett_1_WL_HPa', 'Melian_OO_OO_PSD', 'Dattilo_OO_OO_PSD', 'Dattilo_OO_OO_PA']

corr_importance={}

for name in filtered_network_names:
    my_sign = dict_net_type[dict_name_net[name]]
    Mnet = Read_net_general(name)
    interactions = list(Mnet.slices[2])
    linking_set = get_linking_set(Mnet)
    nodes_to_erase = get_nodes_in_set(Mnet, set_name="Plant")
    N = len(nodes_to_erase)

    filename =  "../../OUTPUT/Data/Node_impact_%s_%s.csv" % (name, "RND")
    Node_impact_df = pd.read_csv(filename, index_col=0)
    Node_impact_df["Robust_merged"]=1-Node_impact_df["Area_merged"]

    if (my_sign == "AA"):  # correct the error in Area parsitism
            # print("CHANGE!###############################################")
            try: #unifying names of interaction
                Node_impact_df["Area_parasitism"] = Node_impact_df["Area_Parasitoid"]
            except:
                Node_impact_df["Area_parasitism"] = Node_impact_df["Area_leaf_miner_parasitoid"]

    rankings_names = ["Area_merged", "Area_%s" % interactions[0], "Area_%s" % interactions[1],
                          "Area_only_%s" % interactions[0], "Area_only_%s" % interactions[1]]
    rankings_names_rev=["Robust_merged","Robust_%s" % interactions[0], "Robust_%s" %interactions[1],"Robust_only_%s" %interactions[0],"Robust_only_%s" % interactions[1]]
    
    for each_interaction in interactions:
        Node_impact_df["Robust_%s" % each_interaction]=1-Node_impact_df["Area_%s" % each_interaction]
        Node_impact_df["Robust_only_%s" % each_interaction]=1-Node_impact_df["Area_only_%s" % each_interaction]
    
    Node_corr_df = Node_impact_df.corr(method="spearman")

    #for each ranking name let's create a concatenating dataframe of correlations
    for rname in rankings_names_rev:
        try:
            corr_importance[rname]
        except:
            corr_importance[rname]=pd.DataFrame() #if doesent exists the first, create it
        corr_importance[rname]=pd.concat([corr_importance[rname],Node_corr_df[rname][:N]])         

types=['Robust_merged', 'Robust_herbivory', 'Robust_pollination','Robust_dispersion','Robust_ant', 'Robust_parasitism' ]


nrow=1
ncol=len(types)
plt.rcParams.update({'font.size': 14})
f, axs = plt.subplots(nrow, ncol, figsize=(5*ncol,4*nrow),sharey=True)
#fontsize = 21      

for i in range(ncol):
    rname=types[i]
    print(i, rname)
    corr_importance[rname][:N].hist(ax=axs[i],bins=50)
    axs[i].set_title(rname)

f.supylabel('N(importance)')
f.supxlabel('importance')
plt.tight_layout()
#plt.show()

#lets create the positive (pollination + dispersion + ant) and negative layers (herbivory + parasitism)
corr_importance['positive']=pd.concat([corr_importance['Robust_pollination'],corr_importance['Robust_dispersion'],corr_importance['Robust_ant']])
corr_importance['negative']=pd.concat([corr_importance['Robust_herbivory'],corr_importance['Robust_parasitism']])

types=['Robust_merged','positive','negative']

rename_dict={"Robust_herbivory":"Herbivore\nranking", "Robust_pollination":"Pollinator\nranking","Robust_merged":"Whole community ranking\n(tri-partite network)",
             "positive": "Positive interactions", "negative":"Negative interactions"}

nrow=1
ncol=3
f, axs = plt.subplots(nrow, ncol, figsize=(5*ncol,4*nrow),sharey=True)
fontsize = 21      

for i in range(ncol):
    rname=types[i]
    print(i, rname)
    corr_importance[rname][:N].hist(ax=axs[i],bins=50)
    axs[i].set_title(rename_dict[rname])

f.supylabel('N(importance)')
f.supxlabel('importance')
plt.tight_layout()
outfilename="../OUTPUT/Images/FIGURE_S17.pdf"
#plt.show()
plt.savefig(outfilename)



