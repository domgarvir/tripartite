from Mnetworks import *
from Graphics import *

ext_MODE="RND"
null_models=[None,"NL", "NL2", "K2" , "K"]
#metric='Area_merged'
metric="Robust_merged"

Nrep=100
ncol=3
nrow=2

f, axarr = plt.subplots(nrow, ncol, figsize=(5*ncol,4*nrow))
axarr=axarr.ravel()
fontsize = 20

index_letter=["A","B","C","x","D","E"]

df_filename = "../OUTPUT/Data/Networks_df_%s.csv" % ext_MODE
emp_df = pd.read_csv(df_filename, index_col=0)
emp_df["Robust_merged"]=1-emp_df["Area_merged"]

index=0
index_text=0
for i in range(len(null_models)):  # null_modesl = row
    if (null_models[i]):
        filename_df_av = "../OUTPUT/Data/Networks_rnd_df_%s_%s_%s_av.csv" % (ext_MODE, null_models[i], Nrep)
        filename_df_std = "../OUTPUT/Data/Networks_rnd_df_%s_%s_%s_std.csv" % (ext_MODE, null_models[i], Nrep)
        rnd_av_df = pd.read_csv(filename_df_av)
        rnd_av_df["Robust_merged"]=1-rnd_av_df["Area_merged"]
        rnd_std_df = pd.read_csv(filename_df_std)
        ylabel = "Zs_%s_%s" % (null_models[i], metric)
        #emp_df.loc[:, ylabel] = calc_Z(emp_df["Area_merged"], rnd_av_df["Area_merged"], emp_df["Area_merged_std"], rnd_std_df["Area_merged"])
        emp_df.loc[:, ylabel] = calc_Z(emp_df["Robust_merged"], rnd_av_df["Robust_merged"], emp_df["Area_merged_std"], rnd_std_df["Area_merged"])
        b = plot_boxplot_fromdf_wseaborn_to_ax(emp_df, ylabel, ax=axarr[index], annot=True, points=True, fontsize=20,ylabel=ylabel)
        axarr[index].axhline(0, color='whitesmoke', linestyle='dotted')
        axarr[index].axhline(1.96, color='gainsboro', linestyle='dotted')
        axarr[index].axhline(-1.96, color='gainsboro', linestyle='dotted')
        axarr[index].axhline(2.33, color='silver', linestyle='dotted')
        axarr[index].axhline(-2.33, color='silver', linestyle='dotted')
        b.tick_params(labelsize=20)
        axarr[index].text(0.01, 0.9, index_letter[index], transform=axarr[index].transAxes, size=20, weight='bold')

        index += 1
        if (index == 3):
            axarr[index].set_axis_off()
            index += 1

    else:
        b = plot_boxplot_fromdf_wseaborn_to_ax(emp_df, metric, ax=axarr[index], annot=True, points=True, ylabel='Robust_merged', fontsize=20, title=ext_MODE)
        axarr[index].text(0.01, 0.9, index_letter[index], transform=axarr[index].transAxes, size=20, weight='bold')
        index += 1
        b.tick_params(labelsize=20)

output_file="../OUTPUT/Images/Figure_S10.pdf"
plt.tight_layout()
plt.savefig(output_file, bbox_inches='tight')