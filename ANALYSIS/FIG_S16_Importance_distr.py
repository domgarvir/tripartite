from Mnetworks import *
from Metrics import *
from Graphics import *

#choose network to plot
name='Sinohara_2_E_PH'#'Sinohara_4_ALL_PH'
ext_MODE="RND"

Mnet=Read_net_general(name)
interactions=list(Mnet.slices[2])
linking_set=get_linking_set(Mnet)
K_df, linking_set = get_full_K_df(name,Mnet,linking_set=linking_set)
LS_df=K_df[K_df["set"]==linking_set]
old_linking_set = get_linking_set(Mnet)

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

#figure
rename_dict={"Area_herbivory":"Herbivore", "Area_pollination":"Pollination","Area_merged":"Whole community\n(tri-partite network)","Robust_herbivory":"Herbivore\nranking", "Robust_pollination":"Pollinator\nranking","Robust_merged":"Whole community ranking\n(tri-partite network)"}

nrow=2 #j; first row with distribution, second with correlation
ncol=3
plt.rcParams.update({'font.size': 18})
f, axs = plt.subplots(nrow, ncol, figsize=(5*ncol,4*nrow),sharey="row", sharex=True)
#fontsize = 21

#we want to do the distribution of importance values

for i in range(ncol):
    rname=rankings_names_rev[i]
    print(i, rname)
    Node_corr_df[rname][:N].hist(ax=axs[0][i],bins=50)
    axs[0][i].set_title(rename_dict[rname])
    axs[0][0].set_ylabel('N(importance)')
   
    #axs[i].hist(Node_corr_df[rname][:N],nbins=20)

    sns.scatterplot(x=Node_corr_df[rname][:N], y=ranking[rname], ax=axs[1][i],  marker="o", s=80, legend=False)
    
    corrfunc(Node_corr_df[rname][:N], -1*ranking[rname], axs[1][i], xy=[0.65,0.4], method="spearman", square=False, fontsize=18)
    corrfunc(Node_corr_df[rname][:N], -1*ranking[rname], axs[1][i], xy=[0.65,0.2], method="pearson", square=False, fontsize=18)

    axs[1][i].set_xlabel('')
    axs[1][0].set_ylabel('Ranking')

plt.gca().invert_yaxis()
#f.supylabel('N(importance)')
f.supxlabel('importance')

#plt.show()
plt.tight_layout()
#f.subplots_adjust(hspace=0.25,wspace=0.2)
outfilename="../OUTPUT/Images/FIGURE_S16.pdf"
#plt.show()
plt.savefig(outfilename)

#correlation with degree 

LS_df["pollination"]
ranking[rname]

nrow=1 #j; first row with distribution, second with correlation
ncol=1
plt.rcParams.update({'font.size': 18})
f, axs = plt.subplots(nrow, ncol, figsize=(6*ncol,4*nrow),sharey="row", sharex=True)

#i=0 #both
rname='Robust_merged'
interaction="O.Degree"
sns.scatterplot(x= LS_df[interaction], y=ranking[rname], ax=axs,  marker="o", s=80, legend=False)
corrfunc(LS_df[interaction], 1*ranking[rname], axs, xy=[0.65,0.4], method="spearman", square=False, fontsize=18)
corrfunc(LS_df[interaction], 1*ranking[rname], axs, xy=[0.65,0.2], method="pearson", square=False, fontsize=18)

axs.set_ylabel('Ranking\n(whole-community)')

plt.gca().invert_yaxis()
plt.tight_layout()
#f.subplots_adjust(hspace=0.25,wspace=0.2)
outfilename="../OUTPUT/Images/FIGURE_R7.pdf"
#plt.show()
plt.savefig(outfilename)