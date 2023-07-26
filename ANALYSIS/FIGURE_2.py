from Mnetworks import *
from Metrics import *
from Graphics import *
pd.options.mode.chained_assignment = None  # default='warn'

#measures to study in the linking set nodes
measures=["cLSsim","CinLHubs_20","new_cLS_PR"]
#measures=["CinLHubs_5","CinLHubs_10","CinLHubs_20"]

#size of figure
#nrow=len(measures)
#ncol=1
ncol=len(measures)
nrow=1
#name of the database dile
df_filename="../OUTPUT/Data/Networks_df_%s.csv" % "RND"
emp_df=pd.read_csv(df_filename)


cm = 1/2.54

f, axarr = plt.subplots(nrow, ncol, figsize=(5.8*ncol,4*nrow))
fontsize = 21

for j in range(ncol):
    #metric=metric_dict[measures[j]]
    plot_boxplot_fromdf_wseaborn_to_ax(emp_df, measures[j], ax=axarr[j], annot=True, points=True, index=string.ascii_uppercase[j], box_by="SIGN", ylabel=measures[j])
    axarr[j].tick_params(labelsize=21)
    #axarr[i][j].set_ylabel(metric, rotation=0, labelpad=15)
    axarr[j].xaxis.label.set_size(21)
    axarr[j].yaxis.label.set_size(19)
    axarr[j].spines['right'].set_visible(False)
    axarr[j].spines['top'].set_visible(False)

axarr[0].set_yticks([0.2,0.4,0.6])
axarr[0].set_ylim([0.0, 0.7])
axarr[1].set_yticks([0.0,0.5,1.0])
axarr[1].set_ylim([0.0, 1.8])
axarr[2].set_yticks([0.4,0.6,0.8,1.0])
axarr[2].set_ylim([0.3, 1.1])

plt.tight_layout()
#f.subplots_adjust(hspace=0.25,wspace=0.2)
outfilename="../OUTPUT/Images/FIGURE_2.pdf"
#plt.show()
plt.savefig(outfilename)