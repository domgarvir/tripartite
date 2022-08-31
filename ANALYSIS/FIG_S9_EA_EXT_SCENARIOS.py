from Mnetworks import *
from Graphics import *

ext_MODES=["RND","OD","ODinv"]

ncol=len(ext_MODES)
nrow=1

f, axarr = plt.subplots(nrow, ncol, figsize=(5*ncol,4*nrow),sharex=True)
fontsize = 20

index=0
#metric='Area_merged'
metric='Robust_merged'
for j in range(len(ext_MODES)): #ext_mode = col
    # load empirical datasets
    df_filename = "../OUTPUT/Data/Networks_df_%s.csv" % ext_MODES[j]
    emp_df = pd.read_csv(df_filename, index_col=0)
    emp_df["Robust_merged"]=1-emp_df["Area_merged"]
    if (j==0):
        ylabel = 'R_merged'
    else:
        ylabel=''
    b = plot_boxplot_fromdf_wseaborn_to_ax(emp_df, metric, ax=axarr[index], annot=True, points=True, ylabel=ylabel, fontsize=20, title=ext_MODES[j])
    axarr[index].text(0.01, 0.9, string.ascii_uppercase[index], transform=axarr[index].transAxes, size=20, weight='bold')
    index +=1
    b.tick_params(labelsize=20)
    b.set_yticks([0.0,0.5,1.0])
    #axarr[index].set_xlabel(fontsize=18)

output_file="../OUTPUT/Images/Figure_S9.pdf"
plt.tight_layout()
plt.savefig(output_file, bbox_inches='tight')