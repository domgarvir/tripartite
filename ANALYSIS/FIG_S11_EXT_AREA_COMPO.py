import semopy as sem
import statsmodels.formula.api as smf
import pandas as pd
#from sklearn import preprocessing
from Graphics import *
from stargazer.stargazer import Stargazer, LineLocation
import imgkit

savename_emp="../OUTPUT/Data/Networks_df_%s.csv" % "RND"
emp_df=pd.read_csv(savename_emp,index_col=0,header=[0])
emp_df["Robust_merged"]=1-emp_df["Area_merged"]
emp_df["R_lA"]=1-emp_df["EA_lA"]
emp_df["R_lB"]=1-emp_df["EA_lB"]
df_M=emp_df[emp_df["sign"] != "AA"]
df_A=emp_df[emp_df["sign"] == "AA"]

####### Now obtaining EA_merged from interaction layer EA
df_M.set_index("name",inplace=True)

#desc1= 'Area_merged ~ EA_lA  +   EA_lB -1'
desc1= 'Robust_merged ~ R_lA  +   R_lB -1'
lm1_M = smf.ols(formula=desc1, data=df_M).fit()
#desc2= 'Area_merged ~ EA_lA -1'
desc2= 'Robust_merged ~ R_lA -1'
lm2_M = smf.ols(formula=desc2, data=df_M).fit()
#desc3= 'Area_merged ~ EA_lB -1'
desc3= 'Robust_merged ~ R_lB -1'
lm3_M = smf.ols(formula=desc3, data=df_M).fit()

#start figure
fig = plt.figure(figsize=(9,7))
gs = fig.add_gridspec(2, 3) #y,x
fontsize = 18

#EA vs EA_larger
ax1 = fig.add_subplot(gs[0:1, 0:1])
#x="EA_lA"
#y="Area_merged"
x="R_lA"
y="Robust_merged"
sns.regplot(x=x, y=y, data=df_M, ax=ax1, scatter=True, marker='v',color='green', line_kws={"linewidth": 1},ci=98)
corrfunc(df_M[x], df_M[y], ax1, method="pearson", color='black', square=True, xy=(0.5,0.2),fontsize=16)
#ax1.set_ylabel("EA merged")
#ax1.set_xlabel("EA large")
ax1.set_ylabel(r"$R$",fontsize=20,rotation=0)
ax1.set_xlabel(r"$R_L$",fontsize=15)
#ax1.set_yticks([0.2,0.3,0.4,0.5])
ax1.set_yticks([0.5,0.6,0.7,0.8])
ax1.text(-0.15, 1.05, "A", transform=ax1.transAxes, size=15, weight='bold')

#EA vs EA_smaller
ax2 = fig.add_subplot(gs[0:1, 1:2])
#x="EA_lB"
x="R_lB"
sns.regplot(x=x, y=y, data=df_M, ax=ax2, scatter=True, marker='^',color='red', line_kws={"linewidth": 1},ci=98)
corrfunc(df_M[x], df_M[y], ax2, method="pearson", color='black', square=True, xy=(0.5,0.2),fontsize=16)
ax2.set_ylabel("")
#ax2.set_yticks([0.2,0.3,0.4,0.5])
ax2.set_yticks([0.5,0.6,0.7,0.8])
#ax2.set_xlabel("EA small")
ax2.set_xlabel(r"$R_S$",fontsize=15)
ax2.text(-0.15, 1.05, "B", transform=ax2.transAxes, size=15, weight='bold')

#EA vs EA_compo
ax4 = fig.add_subplot(gs[0:1, 2:3])
regr_params=lm1_M.params #get regression parameters
#df_M["EA_est"]=regr_params["EA_lA"]*df_M["EA_lA"]+regr_params["EA_lB"]*df_M["EA_lB"]
df_M["R_est"]=regr_params["R_lA"]*df_M["R_lA"]+regr_params["R_lB"]*df_M["R_lB"]
#x="EA_est"
x="R_est"
sns.regplot(x=x, y=y, data=df_M, ax=ax4, scatter=True, marker='o',color='#FA8217', line_kws={"linewidth": 1},ci=98)
corrfunc(df_M[x], df_M[y], ax4, method="pearson", color='black', square=True, xy=(0.5,0.2),fontsize=16)
ax4.set_ylabel("")
#ax4.set_yticks([0.2,0.3,0.4,0.5])
ax4.set_yticks([0.5,0.6,0.7,0.8])
#ax4.set_xlabel("EA est.")
ax4.set_xlabel(r"$R_{(est)}$",fontsize=15)
ax4.text(-0.15, 1.05, "C", transform=ax4.transAxes, size=15, weight='bold')

#EA composition
#df_M.sort_values("Area_merged",inplace=True)
df_M.sort_values("Robust_merged",inplace=True)

ax3 = fig.add_subplot(gs[1:2, 0:3])
#ax3.plot(df_M.index,df_M["Area_merged"],'-o',label="EA merged")
ax3.plot(df_M.index,df_M["Robust_merged"],'-o',label=r"$R$")
##plt.plot(df_M.index,0.604*df_M["EA_lA"]+0.218*df_M["EA_lB"]+0.068,'-o',label="EA stimation")
##plt.plot(df_M.index,0.965*df_M["EA_lA"],'-o',label="EA stimation from EA_lA")
##plt.plot(df_M.index,0.995*df_M["EA_lB"],'-o',label="EA stimation from EA-lB")
#ax3.plot(df_M.index,regr_params["EA_lA"]*df_M["EA_lA"]+regr_params["EA_lB"]*df_M["EA_lB"],'-o',label="EA estimation")
#ax3.plot(df_M.index,df_M["EA_lA"],'v',markersize=3,label="EA large")
#ax3.plot(df_M.index,df_M["EA_lB"],'^',markersize=3,label="EA small")
ax3.plot(df_M.index,regr_params["R_lA"]*df_M["R_lA"]+regr_params["R_lB"]*df_M["R_lB"],'-o',label="R estimation")
ax3.plot(df_M.index,df_M["R_lA"],'v',markersize=3,label="R large")
ax3.plot(df_M.index,df_M["R_lB"],'^',markersize=3,label="R small")
#ax3.annotate(r'$EA_{(est)}$ = {:.2f}'.format(r), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)
#ax3.text(0.4, 0.15, '$EA_{(est)} = %.2f EA_{L} + %.2f EA_{S} $' % (regr_params["EA_lA"], regr_params["EA_lB"]), fontsize=17, transform=ax3.transAxes)
ax3.text(0.6, 0.15, '$R_{(est)} = %.2f R_{L} + %.2f R_{S} $' % (regr_params["R_lA"], regr_params["R_lB"]), fontsize=17, transform=ax3.transAxes)
ax3.tick_params(axis='x', labelrotation= 90)
#ax3.set_ylabel("EA")
ax3.set_ylabel(r"R",rotation=0,fontsize=20)
ax3.legend(loc= 'upper left')
ax3.text(-0.01, 1.05, "D", transform=ax3.transAxes, size=15, weight='bold')

outfilename="../OUTPUT/Images/Figure_S11.pdf"

plt.tight_layout()
plt.savefig(outfilename)


quit()