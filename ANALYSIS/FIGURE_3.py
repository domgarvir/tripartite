from Mnetworks import *
from Metrics import *
from Graphics import *

from sklearn import preprocessing
scaler = preprocessing.StandardScaler()
import statsmodels.formula.api as smf

pd.options.mode.chained_assignment = None  # default='warn'

ext_MODE="RND"

df_RND_filename="../OUTPUT/Data/Networks_df_%s.csv" % ext_MODE
df={}

df["RND"]=pd.read_csv(df_RND_filename, index_col=[0])

#
ncol=3 #len(ext_MODES) # 3 columns
nrow=1 #2 # 3 rows
f, axarr = plt.subplots(nrow, ncol, figsize=(5.8*ncol,4*nrow),sharex=False,sharey=False)
fontsize = 21

str_measures =["Area_merged"]
ylabel=["EA_corr","Area_merged"]
i=0
#correlation extinction area
plot_boxplot_fromdf_wseaborn_to_ax(df[ext_MODE], "r_EA", ax=axarr[i], annot=True, points=True, index=string.ascii_uppercase[i], box_by="SIGN",ylabel=ylabel[i])
axarr[i].tick_params(labelsize=21)
axarr[i].xaxis.label.set_size(21)
axarr[i].yaxis.label.set_size(19)
axarr[i].spines['right'].set_visible(False)
axarr[i].spines['top'].set_visible(False)
axarr[i].set_ylim([-0.30, 1.20])
#axarr[i].set_ylabel[r'Interdependence($I$)']

i=1
df[ext_MODE]["Robust_merged"]= 1 - df[ext_MODE]["Area_merged"]
#plot_boxplot_fromdf_wseaborn_to_ax(df[ext_MODE], "Area_merged", ax=axarr[i], annot=True, index=string.ascii_uppercase[i], points=True, ylabel=ylabel[i], box_by="SIGN")
plot_boxplot_fromdf_wseaborn_to_ax(df[ext_MODE], "Robust_merged", ax=axarr[i], annot=True, index=string.ascii_uppercase[i], points=True, ylabel=ylabel[i], box_by="SIGN")
axarr[i].tick_params(labelsize=21)
axarr[i].xaxis.label.set_size(21)
axarr[i].yaxis.label.set_size(19)
axarr[i].set_yticks([0.0,0.5,1.0])
axarr[i].spines['right'].set_visible(False)
axarr[i].spines['top'].set_visible(False)
#axarr[i].set_ylabel(r"Fragility ($F_M$)")
axarr[i].set_ylabel(r"Robustness ($R$)")

i=2
#include the estimation of the EA
emp_df=df["RND"]

#create columns for robustness as 1-EA:
emp_df["Robust_merged"]=1-emp_df["Area_merged"]
emp_df["R_lA"]=1-emp_df["EA_lA"]
emp_df["R_lB"]=1-emp_df["EA_lB"]

df_M=emp_df[emp_df["sign"] != "AA"] #networks with plants as linking set

#desc1= 'Area_merged ~ EA_lA  +   EA_lB -1'
desc1= 'Robust_merged ~ R_lA  +   R_lB -1'
lm1_M = smf.ols(formula=desc1, data=df_M).fit()
regr_params=lm1_M.params #get regression parameters
#df_M["EA_est"]=regr_params["EA_lA"]*df_M["EA_lA"]+regr_params["EA_lB"]*df_M["EA_lB"]
df_M["Robust_est"]=regr_params["R_lA"]*df_M["R_lA"]+regr_params["R_lB"]*df_M["R_lB"]
#y="EA_est"
x="Robust_est"
#x="Area_merged"
y="Robust_merged"
sns.regplot(x=x, y=y, data=df_M, ax=axarr[i], scatter=False, marker='o',color="#87aa23", line_kws={"linewidth": 2},ci=98)
sns.scatterplot(x=x, y=y, hue="sign", data=df_M, ax=axarr[i],  palette=color_dict, marker="o", s=100, legend=False)
corrfunc(df_M[x], df_M[y], axarr[i], method="pearson", color='black', square=True, xy=(0.2,0.65),fontsize=21)
#axarr[i].text(0.2, 0.15, '$EA_{(est)} = %.2f EA_{L} + %.2f EA_{S} $' % (regr_params["EA_lA"], regr_params["EA_lB"]), fontsize=21, transform=axarr[i].transAxes)
axarr[i].text(0.3, 0.10, '$R_{est} = %.2f R_{L} + %.2f R_{S} $' % (regr_params["R_lA"], regr_params["R_lB"]), fontsize=21, transform=axarr[i].transAxes)
#axarr[i].set_xlabel(r'$EA_{M}$')
#axarr[i].set_ylabel("EA(est)")
#axarr[i].set_ylabel(r"Estimated Fragility ($F_{EST}$)")
#axarr[i].set_xlabel(r"Fragilit ($F_M$)")
axarr[i].set_xlabel(r"Estimated Robustness ($R_{est}$)")
axarr[i].set_ylabel(r"Robustness ($R$)")
#axarr[i].set_yticks([0.2,0.3,0.4,0.5])
axarr[i].text(0.01, 0.9,  "C", transform=axarr[i].transAxes, size=21, weight='black')
axarr[i].tick_params(labelsize=21)
axarr[i].xaxis.label.set_size(19)
axarr[i].yaxis.label.set_size(19)
axarr[i].spines['right'].set_visible(False)
axarr[i].spines['top'].set_visible(False)

outfilename="../OUTPUT/Images/FIGURE_3.pdf"
plt.tight_layout()
#plt.show()
plt.savefig(outfilename)

quit()
