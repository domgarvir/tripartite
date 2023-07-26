import matplotlib.pyplot as plt

from Mnetworks import *
from Metrics import *
from Graphics import *
import matplotlib.gridspec as gridspec
from pathlib import Path

#choose network to plot
name='Sinohara_2_E_PH'#'Sinohara_4_ALL_PH'
ext_MODE="RND"

Mnet=Read_net_general(name)
interactions=list(Mnet.slices[2])
linking_set=get_linking_set(Mnet)
K_df, linking_set = get_full_K_df(name,Mnet,linking_set=linking_set)
plot_gml2(Mnet,K_df,name,color_by="set")
sp_sets=list(Mnet.slices[1])
plants=get_nodes_in_set(Mnet,set_name="Plant")
N=len(plants)

#get impact of each plant
filename="../OUTPUT/Data/Ext_Areas/Node_impact_%s_%s.csv" % (name,ext_MODE)
Node_impact_df = pd.read_csv(filename, index_col=0)
Node_impact_df["Robust_merged"]=1-Node_impact_df["Area_merged"]

#Build rankings dataframe
rankings_names=["Area_merged","Area_%s" % interactions[0], "Area_%s" %interactions[1],"Area_only_%s" %interactions[0],"Area_only_%s" % interactions[1]]
rankings_names_rev=["Robust_merged","Robust_%s" % interactions[0], "Robust_%s" %interactions[1],"Robust_only_%s" %interactions[0],"Robust_only_%s" % interactions[1]]
for each_interaction in interactions:
    Node_impact_df["Robust_%s" % each_interaction]=1-Node_impact_df["Area_%s" % each_interaction]
    Node_impact_df["Robust_only_%s" % each_interaction]=1-Node_impact_df["Area_only_%s" % each_interaction]

ranking={}
#measure coefficient od determination to build ranking
Node_corr_df = Node_impact_df.corr(method="spearman")#only ranking, not searching for lineal dependency

for rname in rankings_names_rev:
    ranking[rname]=Node_corr_df[rname][:N].rank(ascending=False)
    ranking[rname].name=rname
    ranking[rname].to_frame().reset_index()
list_of_rankings=list(ranking.values())
all_ranks=pd.concat(list_of_rankings,axis=1)
#all_ranks=all_ranks.sort_values("Area_merged", ascending=False)
all_ranks=all_ranks.sort_values("Robust_merged", ascending=False)
#choose columns names to compare the different rankings
#columns=["Area_%s" %interactions[0],"Area_merged","Area_%s" %interactions[1]]
columns=["Robust_%s" %interactions[0],"Robust_merged","Robust_%s" %interactions[1]]

#FIGURE:
rename_dict={"Area_herbivory":"Herbivore", "Area_pollination":"Pollination","Area_merged":"Whole community\n(tri-partite network)","Robust_herbivory":"Herbivore\nranking", "Robust_pollination":"Pollinator\nranking","Robust_merged":"Whole community ranking\n(tri-partite network)"}
fig = plt.figure(figsize=(13,6))
gs = fig.add_gridspec(3, 8) #y,x
plt.subplots_adjust(hspace=0.0)
fontsize = 20
#ranking igure
ax1 = fig.add_subplot(gs[:, 2:6],frameon=False)
all_ranks[columns].T.rename(index=rename_dict).plot(legend=False,ax=ax1,cmap='viridis',style='.-',marker='o',markersize=11)
#add labels
df_T=all_ranks[columns].T
species_list=list(df_T)
for sp in species_list:
    x=list(df_T.index)
    y=df_T[sp]
    #addlabels(x,y,[N+1 - df_T[sp]["Area_merged"]]*len(x),fontsize=6.8,fontweight='black')
    addlabels(x,y,[N+1 - df_T[sp]["Robust_merged"]]*len(x),fontsize=6.8,fontweight='black')
ax1.yaxis.set_visible(False)
ax1.tick_params(axis=u'both', which=u'both',length=0)
ax1.text(-0.02, 1, "G", transform=ax1.transAxes, size=12, weight='bold')

ax1.set_title("Plant Ranking")
#correlations # columns= herbivory, merged, pollination
rename_dict_2={"Area_herbivory":r'$EA_{H}$', "Area_pollination":r'$EA_{P}$',"Area_merged":r'$EA_{M}$',"Robust_herbivory":r'$R_{Herb.}$',"Robust_pollination":r'$R_{Pol.}$',"Robust_merged":r'$R$'}
plants_to_plot=["Persicaria hydropiper","Lactuca indica","Cerastium fontanum subsp. vulgare var. angustifolium","Pennisetum alopecuroides"]
###################################################
ax2 = fig.add_subplot(gs[0, 0])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[-1, 0])
axs=[ax2,ax3,ax4]
index=["A","B","C"]
i_plant=0

for layer in range(3): #herb, merged, Poll

    #sns.scatterplot(x=plants_to_plot[i_plant], y=columns[layer],  data=Node_impact_df, ax=axs[layer], marker="o",s=10, color='grey', alpha=0.3)
    sns.regplot(x=plants_to_plot[i_plant], y=columns[layer],  data=Node_impact_df, ax=axs[layer],  scatter_kws={"color": "grey","alpha":0.1,"s":10,"marker":"o"}, line_kws={"color": "black"},  ci=98)

    corrfunc(Node_impact_df[plants_to_plot[i_plant]], Node_impact_df[columns[layer]], axs[layer], method="spearman", xy=[0.1,0.1], square=False, fontsize=13)
    axs[layer].set_ylabel(rename_dict_2[columns[layer]], fontsize=14)
    #axs[layer].set_xlabel("", fontsize=16)
    axs[layer].xaxis.set_visible(False)
    axs[layer].text(0.01, 0.9, index[layer], transform=axs[layer].transAxes, size=14, weight='bold')
ax4.xaxis.set_visible(True)
#ax4.set_xlabel("Ext. order", fontsize=14)
ax4.set_xlabel("Position", fontsize=14)

ax2.set_title("Plant 1")
##################################################
ax5 = fig.add_subplot(gs[0, 1])
ax6 = fig.add_subplot(gs[1, 1])
ax7 = fig.add_subplot(gs[-1, 1])
axs=[ax5,ax6,ax7]
index=["D","E","F"]
i_plant=1
for layer in range(3): #herb, merged, Poll

    #sns.scatterplot(x=plants_to_plot[i_plant], y=columns[layer],  data=Node_impact_df, ax=axs[layer], marker="o",s=10, color='grey', alpha=0.3)
    sns.regplot(x=plants_to_plot[i_plant], y=columns[layer], data=Node_impact_df, ax=axs[layer],scatter_kws={"color": "grey", "alpha": 0.1, "s": 10, "marker": "o"}, line_kws={"color": "black"}, ci=98)
    corrfunc(Node_impact_df[plants_to_plot[i_plant]], Node_impact_df[columns[layer]], axs[layer], method="spearman", xy=[0.1,0.1], square=False, fontsize=13)
    #axs[layer].set_ylabel("", fontsize=14)
    #axs[layer].set_xlabel("", fontsize=16)
    axs[layer].xaxis.set_visible(False)
    axs[layer].yaxis.set_visible(False)
    axs[layer].text(0.01, 0.9, index[layer], transform=axs[layer].transAxes, size=12, weight='bold')
ax7.xaxis.set_visible(True)
#ax7.set_xlabel("Ext. order", fontsize=14)
ax7.set_xlabel("Position", fontsize=14)
ax5.set_title("Plant 2")
##################################################
#Rigth panel:
#try to load the rankings, if does not exists then create them from the Networks_df
filename="../OUTPUT/Data/All_Ranking_%s.csv" % ("RND") #carefull I'm reading from outside
my_file = Path(filename)
if my_file.is_file():
    df_rank_similitude=pd.read_csv(filename,index_col=0)
else:
    df_rank_similitude=get_ranking_df()

#erase autocorrelation
df_clean=df_rank_similitude[df_rank_similitude['level_0']!= df_rank_similitude['level_1']]
#divide by type AA and M
A_nets=df_clean["sign"]=="AA"
M_nets=df_clean["sign"]!="AA"
#retain only the correlation of EA of each int. layer with the EA merged
cross_metrics=((df_clean["level_0"]!="Area_merged") & (df_clean["level_1"]!="Area_merged"))
df=df_clean[~cross_metrics]
#select those defined by only one layer
defined_by_1_layer=df['corr']>0.9
#select the ones not defined by any layer
all_not_defined_by_any_layer=df['corr']<0.7
#build the histogram
nb_of_layers=df[all_not_defined_by_any_layer].groupby("name").count()
networks=list(nb_of_layers[nb_of_layers['level_0']>1].index)
no_defined=df["name"].isin(networks)
#create a column for clasifying into defined by one, none, or mixed
df.loc[:,"Defined"]=0
df.loc[defined_by_1_layer,"Defined"]=-1
df.loc[no_defined,"Defined"]=1
def_df=df.groupby(["sign","name"]).min()
definition_dict={-1:"Determined\nby 1", 0:"Mixed", 1:"Emergent"}
def_df.loc[:,"Defined_name"]=[definition_dict[x] for x in def_df["Defined"]]
def_df.reset_index(inplace=True)
export_filename="../OUTPUT/Data/ranking_clasification.csv"
def_df.to_csv(export_filename)
#plot the histogram
ax8 = fig.add_subplot(gs[:, 6:8])
ax8.yaxis.tick_right()

ax8.yaxis.set_label_position("right")
sns.countplot(x="Defined_name", hue="sign", data=def_df,palette=color_dict,ax=ax8)
ax8.legend(loc='upper right')
ax8.text(0.01, 0.96, "H", transform=ax8.transAxes, size=12, weight='bold')
ax8.set_xlabel("")
ax8.set_ylabel("Count", fontsize=14)
ax8.set_title("Ranking classification", fontsize=14)

filename="../OUTPUT/Images/FIGURE_42.pdf"
plt.tight_layout(pad=0.3)
plt.savefig(filename)

quit()
