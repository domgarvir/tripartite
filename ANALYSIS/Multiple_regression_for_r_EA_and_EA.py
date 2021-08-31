import semopy as sem
import statsmodels.formula.api as smf
import pandas as pd
import seaborn as sns
import graphviz
from sklearn import preprocessing
from Graphics import *
# Create the Scaler object
scaler = preprocessing.StandardScaler()
from stargazer.stargazer import Stargazer, LineLocation
from statsmodels.iolib.summary2 import summary_col
import imgkit
#import re

savename_emp="../OUTPUT/Data/Networks_df_%s.csv" % "RND"
emp_df=pd.read_csv(savename_emp,index_col=0,header=[0])
df_M=emp_df[emp_df["sign"] != "AA"]
df_A=emp_df[emp_df["sign"] == "AA"]


#first standardize values
col_names = emp_df.columns
#erase names of categorical columns
col_names=col_names[8:]
scaled_df_A=scaler.fit_transform(df_A[col_names])
scaled_df_A = pd.DataFrame(scaled_df_A, columns=col_names)

scaled_df_M=scaler.fit_transform(df_M[col_names])
scaled_df_M = pd.DataFrame(scaled_df_M, columns=col_names)

#regression models:
# first multi regression MODEL for the r_EA with five structural features
#desc= 'r_EA ~ LS_HD + rb_k + cLSsim + CinLHubs_20 + new_cLS_PR'
desc= 'r_EA ~ HD + rb_k + cLSsim + CinLHubs_20 + new_cLS_PR'
lm1_AA = smf.ols(formula=desc, data=scaled_df_A).fit()
lm1_AA.summary()

lm1_M = smf.ols(formula=desc, data=scaled_df_M).fit()
lm1_M.summary()

desc_A= 'r_EA ~ cLSsim '
lm2_AA = smf.ols(formula=desc_A, data=scaled_df_A).fit()
desc_M= 'r_EA ~ CinLHubs_20 + new_cLS_PR + rb_k'
lm2_M = smf.ols(formula=desc_M, data=scaled_df_M).fit()


stargazer = Stargazer([lm1_AA,lm1_M,lm2_AA,lm2_M])
stargazer.add_line("AIC",[round(lm1_AA.aic,2),round(lm1_M.aic,2),round(lm2_AA.aic,2),round(lm2_M.aic,2)],LineLocation.FOOTER_TOP)
stargazer.custom_columns(['AA(all)', 'MA \& MM(all)','AA(best)','MA \& MM(best)'], [1, 1, 1, 1])
stargazer.significant_digits(2)
stargazer.show_degrees_of_freedom(False)
stargazer.show_model_numbers(False)
stargazer.dependent_variable_name("$r_{EA}$")
#stargazer.covariate_order(['LS_HD','rb_k','cLSsim','CinLHubs_20','new_cLS_PR'])
#stargazer.covariate_order(['HD','rb_k','cLSsim','CinLHubs_20','new_cLS_PR'])
#stargazer.covariate_order(['HD','cLSsim','CinLHubs_20','new_cLS_PR'])
label_dict={"cLSsim": '$C$',  "CinLHubs_20": '$H_{C}$', "new_cLS_PR":'$PR_{C}$', 'LS_HD':'$\sigma_{k}/<k>_{LS}$','HD':'$\sigma_{k}/<k>$',"rb_k": '$r_{b}$','HD':r'$\sigma_{k}/<k>$'}
stargazer.rename_covariates(label_dict)
stargazer.show_precision=True

table_html=stargazer.render_html()
imgkit.from_string(table_html, '../OUTPUT/Images/Regression_table_r_EA.png')

table_latex=stargazer.render_latex()
table_filename="../OUTPUT/Images/Regression_table_r_EA.tex"
text_file = open(table_filename, "w")
text_file.write(table_latex)
text_file.close()

#Now multiple regression for EA:
#desc= 'Area_merged ~ LS_HD + rb_k + cLSsim + CinLHubs_20 + new_cLS_PR'
desc= 'Area_merged ~ HD + rb_k + cLSsim + CinLHubs_20 + new_cLS_PR'
l2m1_AA = smf.ols(formula=desc, data=scaled_df_A).fit()
l2m1_AA.summary()

l2m1_M = smf.ols(formula=desc, data=scaled_df_M).fit()
l2m1_M.summary()

#desc_A= 'Area_merged ~ LS_HD + cLSsim '
desc_A= 'Area_merged ~ HD + cLSsim'
l2m2_AA = smf.ols(formula=desc_A, data=scaled_df_A).fit()

desc_M= 'Area_merged ~ HD + cLSsim  '
l2m2_M = smf.ols(formula=desc_M, data=scaled_df_M).fit()

stargazer = Stargazer([l2m1_AA,l2m1_M,l2m2_AA,l2m2_M])
stargazer.add_line("AIC",[round(l2m1_AA.aic,2),round(l2m1_M.aic,2),round(l2m2_AA.aic,2),round(l2m2_M.aic,2)],LineLocation.FOOTER_TOP)
stargazer.custom_columns(['AA(all)', 'MA \& MM(all)','AA(best)','MA \& MM(best)'], [1, 1, 1, 1])
stargazer.significant_digits(2)
stargazer.show_degrees_of_freedom(False)
stargazer.show_model_numbers(False)
stargazer.dependent_variable_name("$EA$")
#stargazer.covariate_order(['LS_HD','rb_k','cLSsim','CinLHubs_20','new_cLS_PR'])
#stargazer.covariate_order(['HD','rb_k','cLSsim','CinLHubs_20','new_cLS_PR'])
label_dict={"cLSsim": '$C$',  "CinLHubs_20": '$H_{C}$', "new_cLS_PR":'$PR_{C}$', 'LS_HD':'$\sigma_{k}/<k>_{LS}$','HD':'$\sigma_{k}/<k>$',"rb_k": '$r_{b}$', 'HD':r'$\sigma_{k}/<k>$'}
stargazer.rename_covariates(label_dict)
stargazer.show_precision=True

table_html=stargazer.render_html()
imgkit.from_string(table_html, '../OUTPUT/Images/Regression_table_EA.png')

table_latex=stargazer.render_latex()
table_filename="../OUTPUT/Images/Regression_table_EA.tex"
text_file = open(table_filename, "w")
text_file.write(table_latex)
text_file.close()

#for both EA and r_EA
stargazer = Stargazer([lm2_AA,lm2_M,l2m2_AA,l2m2_M])
stargazer.add_line("AIC",[round(lm2_AA.aic,2),round(lm2_M.aic,2),round(l2m2_AA.aic,2),round(l2m2_M.aic,2)],LineLocation.FOOTER_TOP)
stargazer.custom_columns(['AA(best)', 'MA \& MM(best)','AA(best)','MA \& MM(best)'], [1, 1, 1, 1])
stargazer.significant_digits(2)
stargazer.show_degrees_of_freedom(False)
stargazer.show_model_numbers(False)
#stargazer.show_adj_r2(False)
stargazer.dependent_variable_name("$r_EA$  ;    $EA$")
#stargazer.covariate_order(['LS_HD','cLSsim','CinLHubs_20','new_cLS_PR'])
#stargazer.covariate_order(['HD','rb_k','cLSsim','CinLHubs_20','new_cLS_PR'])
#stargazer.covariate_order(['HD','cLSsim','CinLHubs_20','new_cLS_PR'])
label_dict={"cLSsim": '$C$',  "CinLHubs_20": '$H_{C}$', "new_cLS_PR":'$PR_{C}$', 'LS_HD':'$\sigma_{k}/<k>_{LS}$','HD':'$\sigma_{k}/<k>$',"rb_k": '$r_{b}$'}
stargazer.rename_covariates(label_dict)
stargazer.show_precision=True
table_html=stargazer.render_html()
imgkit.from_string(table_html, '../OUTPUT/Images/Regression_table_both.png')

table_latex=stargazer.render_latex()
table_filename="../OUTPUT/Images/Regression_table_both.tex"
text_file = open(table_filename, "w")
text_file.write(table_latex)
text_file.close()

quit()

ncol=nrow=1

df_M.sort_values("Area_merged",inplace=True)
f, axes =plt.subplots(nrow,ncol, figsize=(7*ncol,5*nrow),sharey=True, sharex=False)
plt.plot(df_M.index,df_M["Area_merged"],'-o',label="EA merged")
#plt.plot(df_M.index,0.604*df_M["EA_lA"]+0.218*df_M["EA_lB"]+0.068,'-o',label="EA stimation")
#plt.plot(df_M.index,0.965*df_M["EA_lA"],'-o',label="EA stimation from EA_lA")
#plt.plot(df_M.index,0.995*df_M["EA_lB"],'-o',label="EA stimation from EA-lB")
plt.plot(df_M.index,0.739*df_M["EA_lA"]+0.246*df_M["EA_lB"],'-o',label="EA stimation")
plt.plot(df_M.index,df_M["EA_lA"],'v',markersize=3,label="EA large")
plt.plot(df_M.index,df_M["EA_lB"],'^',markersize=3,label="EA small")
plt.xticks(rotation=90)
plt.ylabel("EA")
plt.legend(loc= 'lower right')
outfilename="./OUTPUT/Images/Figure_EAmerged_from_EAlayer.pdf"
plt.tight_layout()
#plt.show()
plt.savefig(outfilename)

ncol=2
f, axes =plt.subplots(nrow,ncol, figsize=(7*ncol,5*nrow),sharey=True, sharex=False)
axes[0].scatter(df_M["EA_lA"],df_M["Area_merged"],marker='v')

axes[1].scatter(df_M["EA_lB"],df_M["Area_merged"],marker='^')

axes[0].set_ylabel("EA merged")
axes[0].set_xlabel("EA large")
axes[1].set_xlabel("EA small")
quit()


desc1= 'Area_merged ~ EA_lA  +   EA_lB'
lm1_M = smf.ols(formula=desc1, data=scaled_df_M).fit()
desc2= 'Area_merged ~ EA_lA '
lm2_M = smf.ols(formula=desc2, data=scaled_df_M).fit()
desc3= 'Area_merged ~ EA_lB '
lm3_M = smf.ols(formula=desc3, data=scaled_df_M).fit()
#lm1_M.summary()
stargazer = Stargazer([lm1_M,lm2_M,lm3_M])
stargazer.add_line("AIC",[round(lm1_M.aic,2),round(lm2_M.aic,2),round(lm3_M.aic,2)],LineLocation.FOOTER_TOP)
stargazer.custom_columns(['Two layers', 'Largest layer','Smaller layer'], [1, 1, 1])
table_html=stargazer.render_html()
imgkit.from_string(table_html, './OUTPUT/Images/Regression_table_EA_composition.png')
table_latex=stargazer.render_latex()
table_filename="./OUTPUT/Images/Regression_table_EA_composition.tex"
text_file = open(table_filename, "w")
text_file.write(table_latex)
text_file.close()



plt.plot(scaled_df_M.index, scaled_df_M["Area_merged"],scaled_df_M["EA_lA"],scaled_df_M["EA_lB"],'-o')

plt.plot(scaled_df_M.index,scaled_df_M["Area_merged"],'-o')
plt.plot(scaled_df_M.index,0.7*scaled_df_M["EA_lA"]+0.4*scaled_df_M["EA_lB"],'-o')
plt.plot(scaled_df_M.index,scaled_df_M["EA_lA"],'v',markersize=5)
plt.plot(scaled_df_M.index,scaled_df_M["EA_lB"],'^',markersize=5)




min_area=df_M[["EA_lA","EA_lB"]].min(axis=1)
max_area=df_M[["EA_lA","EA_lB"]].max(axis=1)
df_M["min_area"]=min_area
df_M["max_area"]=max_area
#now wothout normalization
desc1= 'Area_merged ~ max_area  +   min_area'
lm1_M = smf.ols(formula=desc1, data=df_M).fit()
desc2= 'Area_merged ~ max_area '
lm2_M = smf.ols(formula=desc2, data=df_M).fit()
desc3= 'Area_merged ~ min_area '
lm3_M = smf.ols(formula=desc3, data=df_M).fit()
#lm1_M.summary()
stargazer = Stargazer([lm1_M,lm2_M,lm3_M])
stargazer.add_line("AIC",[round(lm1_M.aic,2),round(lm2_M.aic,2),round(lm3_M.aic,2)],LineLocation.FOOTER_TOP)
stargazer.custom_columns(['Two layers', 'Largest layer','Smaller layer'], [1, 1, 1])
table_html=stargazer.render_html()
imgkit.from_string(table_html, './OUTPUT/Images/Regression_table_EA_composition_nonorm_minmax.png')
