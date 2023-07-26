from Mnetworks import *
from Metrics import *
from Graphics import *
import os.path
pd.options.mode.chained_assignment = None  # default='warn'

network_names=['Sinohara_1_ALL_PH', 'Sinohara_2_ALL_PH', 'Sinohara_3_ALL_PH', 'Sinohara_4_ALL_PH', 'Sinohara_ALL_A_PH', 'Sinohara_ALL_E_PH', 'Sinohara_ALL_I_PH', 'Sinohara_2_E_PH', 'Sinohara_3_E_PH', 'Sinohara_4_I_PH', 'Melian_OO_OO_PH', 'Hackett_1_ALL_PH', 'Hackett_2_ALL_PH', 'Hackett_1_S_PH', 'Hackett_1_GL_PH', 'Pocock_OO_OO_PH',  'Melian_OO_OO_HSD',  'McFayden_ALL_A_HPa', 'McFayden_1_A_HPa', 'McFayden_2_A_HPa', 'McFayden_3_A_HPa', 'McFayden_4_A_HPa', 'McFayden_5_A_HPa', 'McFayden_6_A_HPa', 'McFayden_7_A_HPa', 'McFayden_8_A_HPa', 'McFayden_9_A_HPa', 'McFayden_10_A_HPa', 'McFayden_ALL_B_HPa', 'McFayden_1_B_HPa', 'McFayden_2_B_HPa', 'McFayden_3_B_HPa', 'McFayden_4_B_HPa', 'McFayden_5_B_HPa', 'McFayden_6_B_HPa', 'McFayden_7_B_HPa', 'McFayden_8_B_HPa', 'McFayden_9_B_HPa', 'McFayden_10_B_HPa', 'Hackett_1_ALL_HPa', 'Hackett_1_WL_HPa',   'Melian_OO_OO_PSD', 'Dattilo_OO_OO_PSD','Dattilo_OO_OO_PA']

Nrep=200
null_model="NL"

#look if we already have the file with the information of nodes by set
emp_file_path="../OUTPUT/Data/DH_by_sets_df_emp.csv"
rnd_file_path="../OUTPUT/Data/DHZSC_by_sets_df_%s_%s.csv" % (null_model,Nrep)

if (False == os.path.exists(rnd_file_path)):#if does not exist create file
#if (1):
    print("no file, starting simulation")
    df = pd.DataFrame(columns=["name", "type", "sign", "linking_set", "set_name", "DH"])
    df_rnd_zscore = pd.DataFrame(columns=["name", "type", "sign", "linking_set", "set_name", "DH"])

    index = 0
    index_rnd = 0
    for name in network_names:
        print(name)
        Mnet=Read_net_general(name)
        linking_set = get_linking_set(Mnet)
        int_type = dict_name_net[name]
        sign = dict_net_type[int_type]

        K_df, linking_set = get_full_K_df(name, Mnet, linking_set=linking_set)
        sp_sets = K_df["set"].unique()

        emp_sigma=dict.fromkeys(sp_sets, 0)
        for myset in sp_sets:
            selection=K_df["set"]==myset
            emp_sigma[myset] = (K_df[selection].std()/K_df[selection].mean())["O.Degree"]
            df.loc[index] = ([name, int_type, sign, linking_set, myset, emp_sigma[myset]])
            index += 1


        rnd_sigma_av=dict.fromkeys(sp_sets, 0)
        rnd_sigma_av2 = dict.fromkeys(sp_sets, 0)
        for rep in range(Nrep):
            rnd_Mnet = randomize_pymnet(Mnet, method=null_model)
            rnd_linking_set=get_linking_set(rnd_Mnet)
            rnd_K_df, rnd_linking_set = get_full_K_df(name, rnd_Mnet, linking_set=rnd_linking_set)
            for myset in sp_sets:
                selection = rnd_K_df["set"] == myset
                sigma_v = (rnd_K_df[selection].std() / rnd_K_df[selection].mean())["O.Degree"]
                rnd_sigma_av[myset] += sigma_v / Nrep
                rnd_sigma_av2[myset] += pow(sigma_v, 2) / Nrep

        for myset in sp_sets:
            print("rnd set, write to rnd df: %s, inx:%s"  % (myset,index))
            # df_rnd_av.loc[index] = ([name, net_type, sign, linking_set, myset, rnd_sigma_av[myset]])
            # df_rnd_std.loc[index] = ([name, net_type, sign, linking_set, myset, np.sqrt(np.fabs(pow(rnd_sigma_av[myset],2) - rnd_sigma_av2))])
            z_score = (emp_sigma[myset] - rnd_sigma_av[myset]) / (
            np.sqrt(np.fabs(pow(rnd_sigma_av[myset], 2) - rnd_sigma_av2[myset])))
            df_rnd_zscore.loc[index_rnd] = ([name, int_type, sign, linking_set, myset, z_score])
            index_rnd += 1

    #out_filename="./OUTPUT/Data/DH_by_sets_df_emp.csv"
    df.to_csv(emp_file_path)

    #rnd_out_filename="./OUTPUT/Data/DHZSC_by_sets_df_%s_%s.csv" % (null_model,Nrep)
    print("output to %s" % rnd_file_path)
    df_rnd_zscore.to_csv(rnd_file_path)

else:
    print("already done, reading the datafile")
    #out_filename = "./OUTPUT/Data/DH_by_sets_df_emp.csv"
    #rnd_out_filename = "./OUTPUT/Data/DHZSC_by_sets_df_%s_%s.csv" % (null_model, Nrep)
    df=pd.read_csv(emp_file_path,index_col=0)
    df_rnd_zscore=pd.read_csv(rnd_file_path, index_col=0)


#now we should have the dataframes
nrow=2
ncol=1

f, axarr = plt.subplots(nrow, ncol, figsize=(8*ncol,4*nrow),sharex=True)
axarr = axarr.ravel()
fontsize = 20
df.rename(columns={"set_name":"sp. set"},inplace=True)
df_rnd_zscore.rename(columns={"set_name":"sp. set"},inplace=True)
b=plot_boxplot_fromdf_wseaborn_to_ax(df,"DH",box_by= "sp. set",ax=axarr[0],index="A",ylabel= "sigmaK_norm",points=True)
b.tick_params(labelsize=15)
c=plot_boxplot_fromdf_wseaborn_to_ax(df_rnd_zscore,"DH",box_by= "sp. set",index="C",ylabel= "sigmaK_norm_Zscore",points=True, significance_lines=True,ax=axarr[1])
c.tick_params(labelsize=15)
axarr[1].tick_params(labelsize=14,labelrotation=30)

filename="../OUTPUT/Images/Figure_S4A.pdf"
plt.tight_layout()
plt.savefig(filename)



#comparing merged vs set degree heterogeneity - Table S3
emp_df_filename="../OUTPUT/Data/Networks_df_%s.csv" % "RND"
emp_df=pd.read_csv(emp_df_filename, index_col=0)

gta=emp_df["LS_HD"]>emp_df["lsA_HD"]
gtb=emp_df["LS_HD"]>emp_df["lsB_HD"]

A= emp_df[gta | gtb].index #at least LS_HD is greather than one of the layers
B= emp_df[gta & gtb].index # LS_HD is greater in the merged than in any of the components

columns=['name','name_layer_A','name_layer_B',"N_A","N_B","LS_HD","lsA_HD","lsB_HD"]
het_df=emp_df[columns]
het_df["Increaser LS HD"]=''
het_df.loc[B,"Increaser LS HD"]='*'
Table_filename="../OUTPUT/Images/Table_S3.tex"
table_latex=het_df.to_latex(float_format="%.2f")
text_file = open(Table_filename, "w")
text_file.write(table_latex)
text_file.close()


#comparing size of layers
emp_df.loc[:,"DN"]=emp_df["N_A"]- emp_df["N_B"]
emp_df.loc[:,"DN_r"]=(emp_df["N_A"] - emp_df["N_B"])/(emp_df["N_A"])

sns.set_style("white")
nrow=1
ncol=2
f, axarr = plt.subplots(nrow, ncol, figsize=(8*ncol,4*nrow),sharex=True)
axarr = axarr.ravel()
fontsize = 20
b=plot_boxplot_fromdf_wseaborn_to_ax(emp_df,"DN",box_by= "INT",ax=axarr[0], index="A", ylabel= "DN",points=True)
b.tick_params(labelsize=15)

c=plot_boxplot_fromdf_wseaborn_to_ax(emp_df,"DN_r",box_by= "INT",ax=axarr[1], index="B", ylabel= "DN_r",points=True)
c.tick_params(labelsize=15)

filename="../OUTPUT/Images/Figure_R6.pdf"
plt.tight_layout()
plt.savefig(filename)


quit()
