from Mnetworks import *
from Metrics import *
from Graphics import *

ext_MODES=["RND","OD", "ODinv"]


str_metrics2=["HD","rb_k"]
str_metrics3=["cLSsim","CinLHubs_20","new_cLS_PR"]

ncol=len(ext_MODES)

#vs connector nodes metrics
nrow=len(str_metrics3)
f, axarr = plt.subplots(nrow, ncol, figsize=(5*ncol,4*nrow),sharex="row")
fontsize = 20

index=0
for j in range(len(ext_MODES)): #ext_mode = col
    # load empirical datasets
    df_filename = "../OUTPUT/Data/Networks_df_%s.csv" % ext_MODES[j]
    emp_df = pd.read_csv(df_filename, index_col=0)
    emp_df["Robust_merged"]=1-emp_df["Area_merged"]

    if (j==0):
        #ylabel="Area_merged"
        ylabel="Robust_merged"
    else:
        ylabel=""

    if (j==1):
        corr_xy = (0.65, 0.7)
    else:
        corr_xy = (0.05, 0.7)

    for i in range(len(str_metrics3)):
        if (i==0):
            title= ext_MODES[j]
        else:
            title=""
        plot_regplot_fromdf_wseaborn_to_ax(emp_df, x=str_metrics3[i], y="Robust_merged", ax=axarr[i][j],xlabel=str_metrics3[i], ylabel=ylabel, index=string.ascii_uppercase[index],title=title,corr_xy=corr_xy)
        index += 1

outfilename="../OUTPUT/Images/Figure_S13.pdf"
plt.tight_layout()
#plt.show()
plt.savefig(outfilename)

#vs basic structural features
nrow=len(str_metrics2)
f, axarr = plt.subplots(nrow, ncol, figsize=(5*ncol,4*nrow),sharex="row")
fontsize = 20

index=0
corr_xy=corr_xy = (0.65, 0.7)
for j in range(len(ext_MODES)): #ext_mode = col
    # load empirical datasets
    df_filename = "../OUTPUT/Data/Networks_df_%s.csv" % ext_MODES[j]
    emp_df = pd.read_csv(df_filename, index_col=0)
    emp_df["Robust_merged"]=1-emp_df["Area_merged"]

    if (j==0):
        ylabel="Robust_merged"
    else:
        ylabel=""

    for i in range(len(str_metrics2)):
        if (i==0):
            title= ext_MODES[j]
        else:
            title=""
        plot_regplot_fromdf_wseaborn_to_ax(emp_df, x=str_metrics2[i], y="Robust_merged", ax=axarr[i][j],xlabel=str_metrics2[i], ylabel=ylabel, index=string.ascii_uppercase[index],title=title,corr_xy=corr_xy)
        index += 1

outfilename="../OUTPUT/Images/Figure_S12.pdf"
plt.tight_layout()
#plt.show()
plt.savefig(outfilename)