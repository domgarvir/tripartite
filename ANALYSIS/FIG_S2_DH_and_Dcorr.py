from Mnetworks import *
from Metrics import *
from Graphics import *
pd.options.mode.chained_assignment = None  # default='warn'

measures=["HD","r_k"]
null_models=["NL","K"]

Nrep=100

#size of figure
ncol=len(measures)
nrow=2
f, axarr = plt.subplots(nrow, ncol, figsize=(5*ncol,4*nrow),sharex=True)
fontsize = 20

i=0
#first row only the barplot of metrics
emp_df_filename="../OUTPUT/Data/Networks_df_%s.csv" % "RND" #for structural metrics we dont care about the extinction order
emp_df=pd.read_csv(emp_df_filename, index_col="name")

Table=emp_df.loc[:,["sign","int","HD","r_k"]]
index=[["A","B"],["C","D"]]
for j in range(len(measures)):
    b=plot_boxplot_fromdf_wseaborn_to_ax(emp_df, measures[j], ax=axarr[i][j], annot=True, points=True, ylabel=measures[j],fontsize=30)
    b.tick_params(labelsize=15)
    axarr[i][j].text(0.01, 0.9, index[i][j], transform=axarr[i][j].transAxes, size=20, weight='bold')


i=1
#second row the zscore of the metrics
for j in range(2):
    filename_df_av = "../OUTPUT/Data/Networks_rnd_df_%s_%s_%s_av.csv" % ("RND", null_models[j], Nrep)
    filename_df_std = "../OUTPUT/Data/Networks_rnd_df_%s_%s_%s_std.csv" % ("RND", null_models[j], Nrep)

    rnd_av_df = pd.read_csv(filename_df_av, index_col="name")
    rnd_std_df = pd.read_csv(filename_df_std, index_col="name")

    metric = measures[j]
    print(metric)
    print("empirical value:%s" % emp_df.loc[:, metric])
    print("average in null: %s" % rnd_av_df.loc[:, metric])
    # print("std in null: %s" % rnd_std_df.loc[:, metric])
    ylabel = "Zs_%s_%s" % (null_models[j], metric)
    emp_df.loc[:, ylabel] = ((emp_df.loc[:, metric] - rnd_av_df.loc[:, metric]))
    emp_df.loc[:, ylabel] = [0 if np.abs(x) < 0.0000001 else x for x in emp_df.loc[:, ylabel]]
    #rnd_std_df.loc[:, metric] = [0 if np.abs(x) < 0.0000000001 else x for x in rnd_std_df.loc[:, metric]]
    emp_df.loc[:, ylabel] = emp_df.loc[:, ylabel].div(rnd_std_df.loc[:, metric]).replace(np.inf, 0).replace(-np.inf,0).fillna(0)
    Table.loc[:,ylabel]=emp_df.loc[:,ylabel]
    b = plot_boxplot_fromdf_wseaborn_to_ax(emp_df, ylabel, ax=axarr[i][j], annot=True, points=True, fontsize=20,ylabel=ylabel)
    axarr[i][j].axhline(0, color='whitesmoke', linestyle='dotted')
    axarr[i][j].axhline(1.96, color='gainsboro', linestyle='dotted')
    axarr[i][j].axhline(-1.96, color='gainsboro', linestyle='dotted')
    axarr[i][j].axhline(2.33, color='silver', linestyle='dotted')
    axarr[i][j].axhline(-2.33, color='silver', linestyle='dotted')
    b.tick_params(labelsize=15)
    axarr[i][j].text(0.01, 0.9, index[i][j], transform=axarr[i][j].transAxes, size=20, weight='bold')

output_file="../OUTPUT/Images/Figure_S2.pdf"
plt.tight_layout()
plt.savefig(output_file, bbox_inches='tight')

#now change labels and export table :
for net in Table.index:
    for measure in ["Zs_%s_%s" % (null_models[0], measures[0]), "Zs_%s_%s" % (null_models[1],measures[1])]:
        value = Table.loc[net,measure]
        if (np.fabs(value) > 3):
            sig = "**"
        elif (np.fabs(value) > 2):
            sig = "*"
        elif (np.fabs(value) > 1):
            sig = ""
        else:
            sig = ""

        Table.loc[net,measure ]= "%.2f(%s)" % (value, sig)
        Table.loc[net,"Ref."]=cite_dict[net.split("_")[0]]

Table.rename(columns=label_dict, inplace=True)
print(Table)
Table_filename="../OUTPUT/Images/Table_1.tex"
table_latex=Table.to_latex(float_format="%.2f")
text_file = open(Table_filename, "w")
text_file.write(table_latex)
text_file.close()