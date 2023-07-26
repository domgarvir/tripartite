from Mnetworks import *
from Metrics import *
from Graphics import *

import os.path
pd.options.mode.chained_assignment = None  # default='warn'

network_names=['Sinohara_1_ALL_PH', 'Sinohara_2_ALL_PH', 'Sinohara_3_ALL_PH', 'Sinohara_4_ALL_PH', 'Sinohara_ALL_A_PH', 'Sinohara_ALL_E_PH', 'Sinohara_ALL_I_PH', 'Sinohara_2_E_PH', 'Sinohara_3_E_PH', 'Sinohara_4_I_PH', 'Melian_OO_OO_PH', 'Hackett_1_ALL_PH', 'Hackett_2_ALL_PH', 'Hackett_1_S_PH', 'Hackett_1_GL_PH', 'Pocock_OO_OO_PH',  'Melian_OO_OO_HSD',  'McFayden_ALL_A_HPa', 'McFayden_1_A_HPa', 'McFayden_2_A_HPa', 'McFayden_3_A_HPa', 'McFayden_4_A_HPa', 'McFayden_5_A_HPa', 'McFayden_6_A_HPa', 'McFayden_7_A_HPa', 'McFayden_8_A_HPa', 'McFayden_9_A_HPa', 'McFayden_10_A_HPa', 'McFayden_ALL_B_HPa', 'McFayden_1_B_HPa', 'McFayden_2_B_HPa', 'McFayden_3_B_HPa', 'McFayden_4_B_HPa', 'McFayden_5_B_HPa', 'McFayden_6_B_HPa', 'McFayden_7_B_HPa', 'McFayden_8_B_HPa', 'McFayden_9_B_HPa', 'McFayden_10_B_HPa', 'Hackett_1_ALL_HPa', 'Hackett_1_WL_HPa',   'Melian_OO_OO_PSD', 'Dattilo_OO_OO_PSD','Dattilo_OO_OO_PA']

Nrep=200
null_model="K"

#look if we already have the file with the information of nodes by set
emp_file_path="../OUTPUT/Data/r_by_layers_df_emp.csv"
rnd_file_path="../OUTPUT/Data/rZSC_by_layers_df_%s_%s.csv" % (null_model,Nrep)

if (False == os.path.exists(rnd_file_path)):#if does not exist create file

    df = pd.DataFrame(columns=["name", "type", "sign", "linking_set", "set_name", "rb","r"])
    df_rnd_zscore = pd.DataFrame(columns=["name", "type", "sign", "linking_set", "set_name", "rb","r"])

    index = 0
    index_rnd = 0
    for name in network_names:
        print(name)
        Mnet = Read_net_general(name)
        interactions = list(Mnet.slices[2])
        linking_set = get_linking_set(Mnet)
        int_type = dict_name_net[name]
        sign = dict_net_type[int_type]

        r_emp = dict.fromkeys(interactions, 0)
        rb_emp = dict.fromkeys(interactions, 0)
        for interaction in interactions:
            sub_net=get_1Aspect_network(Mnet,interaction)
            r, av_k, av_k2, av_k3 = calc_degree_degree_correlations(sub_net)
            rb, av_kb, av_k2b, av_k3b = calc_degree_degree_correlations_bipart(sub_net,linking_set)
            #print(interaction,r)
            r_emp[interaction]=r
            rb_emp[interaction]=rb
            df.loc[index] = ([name, int_type, sign, linking_set, interaction, rb,r])
            index += 1

        rnd_r_av = dict.fromkeys(interactions, 0)
        rnd_r_av2 = dict.fromkeys(interactions, 0)
        rnd_rb_av = dict.fromkeys(interactions, 0)
        rnd_rb_av2 = dict.fromkeys(interactions, 0)
        for rep in range(Nrep):
            rnd_Mnet = randomize_pymnet(Mnet, method=null_model)
            for interaction in interactions:
                rnd_sub_net = get_1Aspect_network(rnd_Mnet, interaction)
                r, av_k, av_k2, av_k3 = calc_degree_degree_correlations(rnd_sub_net)
                rb, av_kb, av_k2b, av_k3b = calc_degree_degree_correlations_bipart(rnd_sub_net, linking_set)
                rnd_r_av[interaction] += r/Nrep
                rnd_r_av2[interaction] +=  pow(r,2)/Nrep
                rnd_rb_av[interaction] += rb / Nrep
                rnd_rb_av2[interaction] += pow(rb, 2) / Nrep

        for interaction in interactions:
            rnd_std=np.sqrt(np.fabs(rnd_r_av2[interaction] - pow(rnd_r_av[interaction],2)))
            zscore= (r_emp[interaction] - rnd_r_av[interaction]) / rnd_std
            #df_rnd_zscore.loc[index_rnd]=([name, int_type, sign, linking_set, interaction, zscore])
            rnd_stdb = np.sqrt(np.fabs(rnd_rb_av2[interaction] - pow(rnd_rb_av[interaction], 2)))
            zscoreb = (rb_emp[interaction] - rnd_rb_av[interaction]) / rnd_stdb

            df_rnd_zscore.loc[index_rnd] = ([name, int_type, sign, linking_set, interaction,  zscoreb, zscore])
            index_rnd += 1

    df.to_csv(emp_file_path)
    df_rnd_zscore.to_csv(rnd_file_path)

else:
    print("already done, reading the datafile")
    #out_filename = "./OUTPUT/Data/DH_by_sets_df_emp.csv"
    #rnd_out_filename = "./OUTPUT/Data/DHZSC_by_sets_df_%s_%s.csv" % (null_model, Nrep)
    df=pd.read_csv(emp_file_path,index_col=0)
    df_rnd_zscore=pd.read_csv(rnd_file_path, index_col=0)

nrow=2
ncol=1
df.rename(columns={"set_name":"int. layer"},inplace=True)
df_rnd_zscore.rename(columns={"set_name":"int. layer"},inplace=True)
f, axarr = plt.subplots(nrow, ncol, figsize=(8*ncol,4*nrow),sharex=True)
axarr = axarr.ravel()
fontsize = 20
b=plot_boxplot_fromdf_wseaborn_to_ax(df,"r",box_by= "int. layer",ax=axarr[0], index="B", ylabel= "r",points=True,y_lim=[-0.7,0.3])
b.tick_params(labelsize=15)
c=plot_boxplot_fromdf_wseaborn_to_ax(df_rnd_zscore,"r",box_by= "int. layer",ax=axarr[1],index="D", ylabel= "r(k)_Zscore_K", points=True, significance_lines=True)
c.tick_params(labelsize=15)
axarr[1].tick_params(labelsize=15,labelrotation=30)

filename="../OUTPUT/Images/Figure_S4B.pdf"
plt.tight_layout()
plt.savefig(filename)

#quit()



quit()
#now we should have the dataframes
#Empirical
plot_boxplot_fromdf_wseaborn_to_ax(df,"rb",box_by= "set_name",points=True)
plt.show()

#Zscore
plot_boxplot_fromdf_wseaborn_to_ax(df_rnd_zscore,"rb",box_by= "set_name",points=True, significance_lines=True)
plt.show()
