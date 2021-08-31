from Mnetworks import *
from Metrics import *
from Graphics import *


savename_emp="../OUTPUT/Data/Networks_df_%s.csv" % "RND"
emp_df=pd.read_csv(savename_emp,index_col=0,header=[0])

structural_measures=["HD","LS_HD","rb_k","cLSsim","CinLHubs_20","new_cLS_PR"]

nrow=2
ncol=ceil(len(structural_measures)/2)


f, axes =plt.subplots(nrow,ncol, figsize=(5*ncol,4*nrow),sharey=True, sharex=False)
axarr=axes.ravel()


for i in range(len(structural_measures)):
    if not (i%ncol):
        ylabel="EA_corr"
    else:
        ylabel=""

    if (i<3):
        corr_xy=(0.65,0.7)
    else:
        corr_xy = (0.05, 0.7)

    plot_regplot_fromdf_wseaborn_to_ax(emp_df, x=structural_measures[i], y="r_EA", ax=axarr[i], xlabel=structural_measures[i], ylabel= ylabel,index=string.ascii_uppercase[i],corr_xy=corr_xy)

outfilename="../OUTPUT/Images/Figure_S7.pdf"
plt.tight_layout()
#plt.show()
plt.savefig(outfilename)

