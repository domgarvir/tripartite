import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import seaborn as sns
from Graphics import corrfunc

color_dict={'AA': "#8711ac", "MA":  "#87aa23","MM":"#4553c2", "M":"#0099cc"} #this

df_filename="../OUTPUT/Data/Networks_df_%s_W.csv" % "RND"
emp_df=pd.read_csv(df_filename, index_col=0)

#r_EA vs r_EAW
emp_df[["r_EA","r_EAW"]]

# EA vs EAW
emp_df[["Area_merged","AreaW_merged"]]

#reconstruct EAmerged from both EAs
emp_df["Robust_merged"]=1-emp_df["Area_merged"]
emp_df["R_lA"]=1-emp_df["EA_lA"]
emp_df["R_lB"]=1-emp_df["EA_lB"]

emp_df["RobustW_merged"]=1-emp_df["AreaW_merged"]
emp_df["RW_lA"]=1-emp_df["EAW_lA"]
emp_df["RW_lB"]=1-emp_df["EAW_lB"]

df_M=emp_df[emp_df["sign"] != "AA"]
df_A=emp_df[emp_df["sign"] == "AA"]


####### Correlation EA_merged Weighted EA_merged
ncol=2 #len(ext_MODES) # 3 columns
nrow=1 #2 # 3 rows
f, axarr = plt.subplots(nrow, ncol, figsize=(6.2*ncol,4.4*nrow),sharex=False,sharey=False)
fontsize = 19

i=0 #EA
x="Robust_merged"
y="RobustW_merged"
sns.regplot(x=x, y=y, data=df_M, ax=axarr[i], scatter=False, marker='o',color="#87aa23", line_kws={"linewidth": 2},ci=98)
sns.scatterplot(x=x, y=y, hue="sign", data=df_M, ax=axarr[i],  palette=color_dict, marker="o", s=100, legend=False)
corrfunc(df_M[x], df_M[y], axarr[i], method="pearson", color='black', square=True, xy=(0.6,0.2),fontsize=21)
#axarr[i].set_xlim([0.52, 0.62])
#axarr[i].set_ylim([0.45, 0.50])
axarr[i].set_xlabel(r"Robustness ($R$)")
axarr[i].set_ylabel(r"Weighted Robustness ($R'$)")
axarr[i].xaxis.label.set_size(18)
axarr[i].yaxis.label.set_size(18)

i=1 #r_EA
x="r_EA"
y="r_EAW"
sns.regplot(x=x, y=y, data=df_M, ax=axarr[i], scatter=False, marker='o',color="#87aa23", line_kws={"linewidth": 2},ci=98)
sns.scatterplot(x=x, y=y, hue="sign", data=df_M, ax=axarr[i],  palette=color_dict, marker="o", s=100, legend=False)
corrfunc(df_M[x], df_M[y], axarr[i], method="pearson", color='black', square=True, xy=(0.6,0.2),fontsize=21)
axarr[i].set_xlabel(r"Interdependence ($I$)")
axarr[i].set_ylabel(r"Weighted Interdependence ($I'$)")
axarr[i].xaxis.label.set_size(17)
axarr[i].yaxis.label.set_size(17)

outfilename="../OUTPUT/Images/FIGURE_S14.png"
plt.tight_layout()
plt.savefig(outfilename,bbox_inches="tight")


####### Now obtaining EA_merged from interaction layer EA
df_M.set_index("name",inplace=True)

ncol=2 #len(ext_MODES) # 3 columns
nrow=1 #2 # 3 rows
f, axarr = plt.subplots(nrow, ncol, figsize=(6.2*ncol,4.4*nrow),sharex=False,sharey=False)
fontsize = 19

i=0
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

i=1
desc2= 'RobustW_merged ~ RW_lA  +   RW_lB -1'
lm2_M = smf.ols(formula=desc2, data=df_M).fit()
regr_paramsW=lm2_M.params #get regression parameters
##df_M["EA_est"]=regr_params["EA_lA"]*df_M["EA_lA"]+regr_params["EA_lB"]*df_M["EA_lB"]
df_M["RobustW_est"]=regr_paramsW["RW_lA"]*df_M["RW_lA"]+regr_paramsW["RW_lB"]*df_M["RW_lB"]
#y="EA_est"
x="RobustW_est"
#x="Area_merged"
y="RobustW_merged"
sns.regplot(x=x, y=y, data=df_M, ax=axarr[i], scatter=False, marker='o',color="#87aa23", line_kws={"linewidth": 2},ci=98)
sns.scatterplot(x=x, y=y, hue="sign", data=df_M, ax=axarr[i],  palette=color_dict, marker="o", s=100, legend=False)
corrfunc(df_M[x], df_M[y], axarr[i], method="pearson", color='black', square=True, xy=(0.2,0.65),fontsize=21)
axarr[i].set_xlabel(r"Estimated Weighted Robustness ($R'_{est}$)")
axarr[i].set_ylabel(r"Weighted Robustness ($R'$)")
##axarr[i].set_yticks([0.2,0.3,0.4,0.5])
axarr[i].text(0.01, 0.9,  "C", transform=axarr[i].transAxes, size=21, weight='black')
axarr[i].tick_params(labelsize=21)
axarr[i].xaxis.label.set_size(19)
axarr[i].yaxis.label.set_size(19)
axarr[i].spines['right'].set_visible(False)
axarr[i].spines['top'].set_visible(False)

outfilename="../OUTPUT/Images/FIGURE_S15.png"
plt.tight_layout()
plt.savefig(outfilename,bbox_inches="tight")
#plt.show()
