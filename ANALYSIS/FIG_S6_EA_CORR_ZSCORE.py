from Mnetworks import *
from Graphics import *

ext_MODE="RND"
null_models=[None,"NL", "NL2", "K2" , "K"]
metric="r_EA"

Nrep=100
ncol=3
nrow=2

f, axarr = plt.subplots(nrow, ncol, figsize=(5*ncol,4*nrow))
axarr=axarr.ravel()
fontsize = 20

index_letter=["A","B","C","x","D","E"]

#load empiric database:
df_filename="../OUTPUT/Data/Networks_df_%s.csv" % ext_MODE
emp_df=pd.read_csv(df_filename,index_col=0)

index=0
index_text=0
for i in range(len(null_models)):  # null_modesl = row
    if (null_models[i]):
        #load database of randomizations
        filename_df_av = "../OUTPUT/Data/Networks_rnd_df_%s_%s_%s_av.csv" % (ext_MODE, null_models[i], Nrep)
        filename_df_std = "../OUTPUT/Data/Networks_rnd_df_%s_%s_%s_std.csv" % (ext_MODE, null_models[i], Nrep)
        rnd_av_df=pd.read_csv(filename_df_av)
        rnd_std_df = pd.read_csv(filename_df_std)
        ylabel = "Zs_%s_%s" % (null_models[i], metric)
        emp_df.loc[:, ylabel] = ((emp_df.loc[:, metric] - rnd_av_df.loc[:, metric]))
        emp_df.loc[:, ylabel] = [0 if np.abs(x) < 0.0000001 else x for x in emp_df.loc[:, ylabel]]
        rnd_std_df.loc[:, metric] = [0 if np.abs(x) < 0.0000000001 else x for x in rnd_std_df.loc[:, metric]]
        emp_df.loc[:, ylabel] = emp_df.loc[:, ylabel].div(rnd_std_df.loc[:, metric]).replace(np.inf, 0).replace(-np.inf,0).fillna(0)
        b = plot_boxplot_fromdf_wseaborn_to_ax(emp_df, ylabel, ax=axarr[index], annot=True, points=True, ylabel=ylabel,fontsize=20)
        axarr[index].axhline(0, color='whitesmoke', linestyle='dotted')
        axarr[index].axhline(1.96, color='gainsboro', linestyle='dotted')
        axarr[index].axhline(-1.96, color='gainsboro', linestyle='dotted')
        axarr[index].axhline(2.33, color='silver', linestyle='dotted')
        axarr[index].axhline(-2.33, color='silver', linestyle='dotted')
        b.tick_params(labelsize=20)
        axarr[index].text(0.01, 0.9, index_letter[index], transform=axarr[index].transAxes, size=20, weight='bold')

        index +=1
        if (index == 3):
            axarr[index].set_axis_off()
            index += 1

    else:
        b=plot_boxplot_fromdf_wseaborn_to_ax(emp_df, metric, ax=axarr[index], annot=True, points=True, ylabel="EA_corr", index=string.ascii_uppercase[index_text],fontsize=20)
        b.tick_params(labelsize=20)
        axarr[index].text(0.01, 0.9, index_letter[index], transform=axarr[index].transAxes, size=20, weight='bold')
        index += 1


output_file="../OUTPUT/Images/Figure_S6.pdf"
plt.tight_layout()
plt.savefig(output_file, bbox_inches='tight')

