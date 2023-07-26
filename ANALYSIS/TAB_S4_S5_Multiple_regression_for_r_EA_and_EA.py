#import semopy as sem
import statsmodels.formula.api as smf
import pandas as pd
import seaborn as sns
#import graphviz
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
emp_df["Robust_merged"]=1-emp_df["Area_merged"]
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

##### new code
variables=['HD' , 'cLSsim' , 'CinLHubs_20', 'new_cLS_PR'] #, 'rb_k' łhy left this out?????
from itertools import combinations
c1 = [['HD'] , ['cLSsim'] , ['CinLHubs_20'], ['new_cLS_PR']] #, ['rb_k']
c2 = [i for i in combinations(variables,2)]
c3 = [i for i in combinations(variables,3)]
c4 = [i for i in combinations(variables,4)]
#c5 = [i for i in combinations(variables,5)]

my_combinations=[c1,c2,c3,c4]#,c5
min_AIC_I_AA=100
min_AIC_I_M=100
min_AIC_R_AA=100
min_AIC_R_M=100

#generate the different models
for index in range(len(my_combinations)): #c1,c2...
    print("combinations with %s elements" % index)
    for element in my_combinations[index]:
        if (len(element) <2 ): #if only one variable
            var_desc=element[0]
        else:
            var_desc=element[0]
            for var_index in range(1,len(element)):
                var_desc = var_desc + " + " + element[var_index]

        desc_I="r_EA ~ %s" % var_desc
        desc_R="Robust_merged ~ %s" % var_desc

        print(desc_I)
        print(desc_R)

        lm1_AA = smf.ols(formula=desc_I, data=scaled_df_A).fit()
        lm1_AA_AIC=lm1_AA.aic
        lm1_M = smf.ols(formula=desc_I, data=scaled_df_M).fit()
        lm1_M_AIC=lm1_M.aic

        lm2_AA = smf.ols(formula=desc_R, data=scaled_df_A).fit()
        lm2_AA_AIC=lm2_AA.aic
        lm2_M = smf.ols(formula=desc_R, data=scaled_df_M).fit()
        lm2_M_AIC=lm2_M.aic

        #find smaller AIC
        if (lm1_AA_AIC < min_AIC_I_AA):
            min_AIC_I_AA = lm1_AA_AIC
            var_I_AA= desc_I
            best_lm1_AA=lm1_AA
            print ("found new min %s rEA_AA: %s" % (min_AIC_I_AA,desc_I))
        if (lm1_M_AIC < min_AIC_I_M):
            min_AIC_I_M = lm1_M_AIC
            var_I_M = desc_I
            best_lm1_M=lm1_M
            print("found new min rEA_M: %s" % desc_I)
        if (lm2_AA_AIC < min_AIC_R_AA):
            min_AIC_R_AA = lm2_AA_AIC
            var_R_AA= desc_R
            best_lm2_AA=lm2_AA
            print("found new min EA_AA: %s" % desc_R)
        if (lm2_M_AIC < min_AIC_R_M):
            min_AIC_R_M = lm2_M_AIC
            var_R_M = desc_R
            best_lm2_M=lm2_M
            print("found new min EA_M: %s" % desc_R)

print("The minimum AIC are for:")
print("I_AA: %s // %s " % (min_AIC_I_AA,var_I_AA))
print("I_M: %s // %s " % (min_AIC_I_M,var_I_M))
print("R_AA: %s  // %s " % (min_AIC_R_AA,var_R_AA))
print("R_M: %s // %s" % (min_AIC_R_M,var_R_M))

#Now print table
stargazer = Stargazer([best_lm1_AA,best_lm1_M,best_lm2_AA,best_lm2_M])
stargazer.add_line("AIC",[round(best_lm1_AA.aic,2),round(best_lm1_M.aic,2),round(best_lm2_AA.aic,2),round(best_lm2_M.aic,2)],LineLocation.FOOTER_TOP)
stargazer.custom_columns(['AA(best)', 'MA \& MM(best)','AA(best)','MA \& MM(best)'], [1, 1, 1, 1])
stargazer.significant_digits(2)
stargazer.show_degrees_of_freedom(False)
stargazer.show_model_numbers(False)
label_dict={"cLSsim": '$C$',  "CinLHubs_20": '$H_{C}$', "new_cLS_PR":'$PR_{C}$', 'LS_HD':'$\sigma_{k}/<k>_{LS}$','HD':'$\sigma_{k}/<k>$',"rb_k": '$r_{b}$','HD':r'$\sigma_{k}/<k>$'}
stargazer.covariate_order(['HD','cLSsim','CinLHubs_20','new_cLS_PR']) #,'rb_k'
stargazer.rename_covariates(label_dict)
stargazer.show_precision=True
stargazer.dependent_variable_name("$I$  ;    $R$")

table_html=stargazer.render_html()
imgkit.from_string(table_html, '../OUTPUT/Images/Table_1.png')

table_latex=stargazer.render_latex()
table_filename="../OUTPUT/Images/Table_1.tex"
text_file = open(table_filename, "w")
text_file.write(table_latex)
text_file.close()

#getting supp tables: 
#interdependence
desc_I_all= 'r_EA ~ HD +  cLSsim + CinLHubs_20 + new_cLS_PR' #regression model of interdependence on all strutural metrics + rb_k
lm_I_AA = smf.ols(formula=desc_I_all, data=scaled_df_A).fit() #fit for Antagonistic-antagonistic
lm_I_M  = smf.ols(formula=desc_I_all, data=scaled_df_M).fit() #fit for mutualistic-mutualistic and mutualistic-antagonistic
#export table
stargazer = Stargazer([lm_I_AA,lm_I_M,best_lm1_AA,best_lm1_M])
stargazer.add_line("AIC",[round(lm_I_AA.aic,2),round(lm_I_M.aic,2),round(best_lm1_AA.aic,2),round(best_lm1_M.aic,2)],LineLocation.FOOTER_TOP)
stargazer.custom_columns(['AA(all)', 'MA \& MM(all)','AA(best)','MA \& MM(best)'], [1, 1, 1, 1])
stargazer.significant_digits(2)
stargazer.show_degrees_of_freedom(False)
stargazer.show_model_numbers(False)
stargazer.dependent_variable_name("$I$ ")
stargazer.covariate_order(['HD','cLSsim','CinLHubs_20','new_cLS_PR']) #,'rb_k'
stargazer.rename_covariates(label_dict)

table_html=stargazer.render_html()
imgkit.from_string(table_html, '../OUTPUT/Images/Table_S4.png')

table_latex=stargazer.render_latex()
table_filename="../OUTPUT/Images/Table_S4.tex"
text_file = open(table_filename, "w")
text_file.write(table_latex)
text_file.close()

#robustness
desc_R_all= 'Robust_merged ~ HD  + cLSsim + CinLHubs_20 + new_cLS_PR' #regression model of robustness on all structural metrics + rb_k
lm_R_AA = smf.ols(formula=desc_R_all, data=scaled_df_A).fit() #fit for Antagonistic-antagonistic
lm_R_M  = smf.ols(formula=desc_R_all, data=scaled_df_M).fit() #fit for mutualistic-mutualistic and mutualistic-antagonistic
#export table
stargazer = Stargazer([lm_R_AA,lm_R_M,best_lm2_AA,best_lm2_M])
stargazer.add_line("AIC",[round(lm_R_AA.aic,2),round(lm_R_M.aic,2),round(best_lm2_AA.aic,2),round(best_lm2_M.aic,2)],LineLocation.FOOTER_TOP)
stargazer.custom_columns(['AA(all)', 'MA \& MM(all)','AA(best)','MA \& MM(best)'], [1, 1, 1, 1])
stargazer.significant_digits(2)
stargazer.show_degrees_of_freedom(False)
stargazer.show_model_numbers(False)
stargazer.dependent_variable_name("R")
stargazer.covariate_order(['HD','cLSsim','CinLHubs_20','new_cLS_PR']) #,'rb_k'
stargazer.rename_covariates(label_dict)

table_html=stargazer.render_html()
imgkit.from_string(table_html, '../OUTPUT/Images/Table_S5.png')

table_latex=stargazer.render_latex()
table_filename="../OUTPUT/Images/Table_S5.tex"
text_file = open(table_filename, "w")
text_file.write(table_latex)
text_file.close()

##### end new code
quit()

