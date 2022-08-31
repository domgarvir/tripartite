from Mnetworks import *
from Metrics import *
from Graphics import *


null_models=[None,"NL","NL2", "K2","K"]
null_models_names=[None,"NL","NL2","K2","K"]
measures=["HD","LS_HD","r_k","cLSsim","CinLHubs_20","new_cLS_PR"]#

Nrep=100

ncol=len(measures)
nrow=len(null_models)

f, axarr = plt.subplots(nrow, ncol, figsize=(5*ncol,3.8*nrow),sharex=True)
fontsize = 10

emp_df_filename="../OUTPUT/Data/Networks_df_%s.csv" % "RND" #for structural metrics we dont care about the extinction order
emp_df=pd.read_csv(emp_df_filename, index_col="name")


index=0
for i in range(len(null_models)): # 1 null model per row, starting with the empirical values
    model=null_models[i]
    print(model)

    if (model):
        filename_df_av= "../OUTPUT/Data/Networks_rnd_df_%s_%s_%s_av.csv"  % ("RND",null_models[i],Nrep)
        print(filename_df_av)
        filename_df_std= "../OUTPUT/Data/Networks_rnd_df_%s_%s_%s_std.csv" % ("RND",null_models[i],Nrep)
        print(filename_df_std)
        rnd_av_df=pd.read_csv(filename_df_av,index_col="name")
        rnd_std_df=pd.read_csv(filename_df_std,index_col="name")

        for j in range(len(measures)):
            metric=measures[j]
            print(metric)
            print("empirical value:%s" % emp_df.loc[:,metric])
            print("average in null: %s" % rnd_av_df.loc[:,metric])
            #print("std in null: %s" % rnd_std_df.loc[:, metric])
            ylabel="Zs_%s_%s" %(null_models[i],metric)
            emp_df.loc[:,ylabel]=((emp_df.loc[:,metric]-rnd_av_df.loc[:,metric]))
            emp_df.loc[:,ylabel]=[0 if np.abs(x) <0.0000001 else x for x in emp_df.loc[:,ylabel]]
            rnd_std_df.loc[:, metric] = [0 if np.abs(x) < 0.0000000001 else x for x in rnd_std_df.loc[:, metric]]
            emp_df.loc[:,ylabel]=emp_df.loc[:,ylabel].div(rnd_std_df.loc[:,metric]).replace(np.inf, 0).replace(-np.inf, 0).fillna(0)
            b=plot_boxplot_fromdf_wseaborn_to_ax(emp_df, ylabel, ax=axarr[i][j], annot=True, points=True,fontsize=20)
            axarr[i][j].axhline(0, color='whitesmoke', linestyle='dotted')
            axarr[i][j].axhline(1.96, color='gainsboro', linestyle='dotted')
            axarr[i][j].axhline(-1.96, color='gainsboro', linestyle='dotted')
            axarr[i][j].axhline(2.33, color='silver', linestyle='dotted')
            axarr[i][j].axhline(-2.33, color='silver', linestyle='dotted')
            b.tick_params(labelsize=20)

    else: #no null model then only plot empirical boxplots
        for j in range(len(measures)):
            b=plot_boxplot_fromdf_wseaborn_to_ax(emp_df, measures[j], ax=axarr[i][j], annot=True, points=True, title=measures[j],fontsize=30)
            b.tick_params(labelsize=25)
            #index += 1

    # include index by row
    axarr[i][0].text(-0.15, 1.05, string.ascii_uppercase[i], transform=axarr[i][0].transAxes, size=18, weight='bold')
    if(i!=0):
        axarr[i][0].set_ylabel(null_models_names[i])
    index += 1


output_file="../OUTPUT/Images/Figure_S6.pdf"
plt.tight_layout()
plt.savefig(output_file, bbox_inches='tight')

