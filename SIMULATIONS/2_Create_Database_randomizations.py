from Mnetworks import *
from Metrics import *
from Robustness import *
import sys
pd.options.mode.chained_assignment = None  # default='warn'
from pathlib import Path

network_names=['Sinohara_1_ALL_PH', 'Sinohara_2_ALL_PH', 'Sinohara_3_ALL_PH', 'Sinohara_4_ALL_PH', 'Sinohara_ALL_A_PH', 'Sinohara_ALL_E_PH', 'Sinohara_ALL_I_PH', 'Sinohara_2_E_PH', 'Sinohara_3_E_PH', 'Sinohara_4_I_PH', 'Melian_OO_OO_PH', 'Hackett_1_ALL_PH', 'Hackett_2_ALL_PH', 'Hackett_1_S_PH', 'Hackett_1_GL_PH', 'Pocock_OO_OO_PH',  'Melian_OO_OO_HSD',  'McFayden_ALL_A_HPa', 'McFayden_1_A_HPa', 'McFayden_2_A_HPa', 'McFayden_3_A_HPa', 'McFayden_4_A_HPa', 'McFayden_5_A_HPa', 'McFayden_6_A_HPa', 'McFayden_7_A_HPa', 'McFayden_8_A_HPa', 'McFayden_9_A_HPa', 'McFayden_10_A_HPa', 'McFayden_ALL_B_HPa', 'McFayden_1_B_HPa', 'McFayden_2_B_HPa', 'McFayden_3_B_HPa', 'McFayden_4_B_HPa', 'McFayden_5_B_HPa', 'McFayden_6_B_HPa', 'McFayden_7_B_HPa', 'McFayden_8_B_HPa', 'McFayden_9_B_HPa', 'McFayden_10_B_HPa', 'Hackett_1_ALL_HPa', 'Hackett_1_WL_HPa',   'Melian_OO_OO_PSD','Dattilo_OO_OO_PSD','Dattilo_OO_OO_PA']#,

network_name_type_tuples=[]
for name in network_names:
     network_name_type_tuples.append((dict_net_type[dict_name_net[name]], dict_name_net[name], name))
#create index
my_index = pd.MultiIndex.from_tuples(network_name_type_tuples,names=['sign','int','name'])

#define the extinction scenario we want to study

#determine the structural metrics we want to include in the database
measures=["name_set_a","name_set_b","linking_set","name_layer_A","name_layer_B","Na","Nb","Nls","N_lsA","N_lsB","N_A","N_B","K_A","K_B","C_A","C_B", "HD","LS_HD","nonLS_HD", "P_HD", "nonP_HD","lsA_HD","lsB_HD", "r_k","rb_k", "cLSsim","CinLHubs_10","CinLHubs_20","CinLHubs_5","new_cLS_PR","PRofCinLhubs", "Area_merged", "EA_a","EA_b","EA_lA","EA_lB","r_EA"]

all_nets_av_data=pd.DataFrame(index=my_index, columns=measures)
all_nets_std_data=pd.DataFrame(index=my_index, columns=measures)

Nrep=100
Nseq=1#we are not using this to obtain EA simultaneously for now

#ext_MODE="RND"
#null_model="K2"  #(NL), (NL2), (K2), (K)
null_model=sys.argv[2]
ext_MODE=sys.argv[1]


for sign, int, name in network_name_type_tuples:
    print(name)
    Mnet = Read_net_general(name)
    linking_set = get_linking_set(Mnet)
    K_df, linking_set = get_full_K_df(name, Mnet, linking_set=linking_set)
    LS_df = K_df[K_df["set"] == linking_set]
    old_linking_set = get_linking_set(Mnet)

    #########################################################
    # get size of subnetworks
    sub_net_sizes = {}
    sub_net_links = {}
    max_size = 0
    max_name = 0
    interactions = list(Mnet.slices[2])
    for interaction in interactions:
        nodes, links = get_nodes_in_layer(Mnet, interaction_layer=interaction, nlinks=True)
        size = len(nodes)
        sub_net_sizes[interaction] = size
        sub_net_links[interaction] = links
        if (size > max_size):
            max_size = size
            max_name = interaction

    name_layer_A = max_name
    N_A = sub_net_sizes.pop(name_layer_A)
    L_A = sub_net_links.pop(name_layer_A)
    K_A = L_A / N_A

    name_layer_B = list(sub_net_sizes.keys())[0]
    N_B = sub_net_sizes.pop(name_layer_B)
    L_B = sub_net_links.pop(name_layer_B)
    K_B = L_B / N_B

    # get sizes of sets
    sp_sets = list(K_df["set"].unique())

    # get sizes os linking set species in different layers
    #print("sign : %s" % sign)

    nodes_lsA = get_nodes_in_spset(Mnet, sign, species_set=linking_set, interaction_layer=name_layer_A)
    N_lsA = len(nodes_lsA)
    nodes_lsB = get_nodes_in_spset(Mnet, sign, species_set=linking_set, interaction_layer=name_layer_B)
    N_lsB = len(nodes_lsB)
    Nls = len(get_nodes_in_spset(Mnet, sign, species_set=linking_set))
    #print(linking_set, name_layer_A,name_layer_B)

    species_set_a = dict_sp_set_name_from_interaction[sign][name_layer_A]
    species_set_b = dict_sp_set_name_from_interaction[sign][name_layer_B]
    Na = len(get_nodes_in_spset(Mnet, sign, species_set=species_set_a, interaction_layer=name_layer_A))
    Nb = len(get_nodes_in_spset(Mnet, sign, species_set=species_set_b, interaction_layer=name_layer_B))

    #averages for metrics
    filename_df = "../OUTPUT/Data/NEWSEM_df_%s_%s_%s_%s_noEA.csv" % (name, "RND", null_model, Nrep)#work with the RND since the structural part is the same anyway
    my_file = Path(filename_df)

    if my_file.is_file():#if exists just read from the file
        df=pd.read_csv(filename_df,index_col=0)
        mean=df.mean()
        std = df.std()
        colnames=list(mean.index)
        other_colnames=['name_set_a',  'name_set_b',  'linking_set',  'name_layer_A',  'name_layer_B']
        all_nets_av_data.loc[(sign, int, name),colnames]=mean
        all_nets_av_data.loc[(sign, int, name),other_colnames]=df.loc[1,other_colnames]
        all_nets_std_data.loc[(sign, int, name),colnames]=std
        all_nets_std_data.loc[(sign, int, name),other_colnames]=df.loc[1,other_colnames]
    else:
        #sizes and connectivities
        df = pd.DataFrame(index=np.arange(1, Nrep), columns=measures)
        for rep in range(1,Nrep):
            print("rep:%s" % rep)
            rnd_Mnet = randomize_pymnet(Mnet, method=null_model)
            rnd_linking_set = get_linking_set(rnd_Mnet)
            rnd_K_df, rnd_linking_set = get_full_K_df(name,rnd_Mnet,linking_set=rnd_linking_set)
            rnd_LS_df = rnd_K_df[rnd_K_df["set"] == rnd_linking_set]
            old_linking_set = get_linking_set(rnd_Mnet)

            ##################################################
            # sizes and connectivities
            sub_net_sizes = {}
            sub_net_links = {}
            interactions = list(rnd_Mnet.slices[2])
            for interaction in interactions:
                nodes, links = get_nodes_in_layer(rnd_Mnet, interaction_layer=interaction, nlinks=True)
                size = len(nodes)
                sub_net_sizes[interaction] = size
                sub_net_links[interaction] = links

            N_A = sub_net_sizes.pop(name_layer_A)
            L_A = sub_net_links.pop(name_layer_A)
            K_A = L_A / N_A

            N_B = sub_net_sizes.pop(name_layer_B)
            L_B = sub_net_links.pop(name_layer_B)
            K_B = L_B / N_B

            nodes_lsA = get_nodes_in_spset(rnd_Mnet, sign, species_set=rnd_linking_set, interaction_layer=name_layer_A)
            N_lsA = len(nodes_lsA)
            nodes_lsB = get_nodes_in_spset(rnd_Mnet, sign, species_set=rnd_linking_set, interaction_layer=name_layer_B)
            N_lsB = len(nodes_lsB)

            Nls = len(get_nodes_in_spset(rnd_Mnet, sign, species_set=rnd_linking_set))
            Na = len(get_nodes_in_spset(rnd_Mnet, sign, species_set=species_set_a, interaction_layer=name_layer_A))
            Nb = len(get_nodes_in_spset(rnd_Mnet, sign, species_set=species_set_b, interaction_layer=name_layer_B))

            C_A = (L_A * 0.5) / (Na * N_lsA)
            C_B = (L_B * 0.5) / (Nb * N_lsB)

            # Degree heterogeneities and correlations
            ###################################################################
            HD = calc_DH(rnd_K_df)
            #print("rn HD:%s " % (HD))
            LS_HD = calc_DH(rnd_LS_df)
            nonLS_HD = calc_DH(rnd_K_df[rnd_K_df["set"] != linking_set])

            P_HD = calc_DH(rnd_K_df[rnd_K_df["set"] == "Plant"])
            nonP_HD = calc_DH(rnd_K_df[rnd_K_df["set"] != "Plant"])/Nrep

            ls_A_HD = (rnd_K_df.loc[nodes_lsA].std() / rnd_K_df.loc[nodes_lsA].mean())[name_layer_A]
            ls_B_HD = (rnd_K_df.loc[nodes_lsB].std() / rnd_K_df.loc[nodes_lsB].mean())[name_layer_B]

            rb = calc_measure(rnd_Mnet, old_linking_set, rnd_LS_df, "rb_k", K_df=rnd_K_df)
            r = calc_measure(rnd_Mnet,old_linking_set,rnd_LS_df,"r_k",K_df=rnd_K_df)

            # Linking set metrics
            #################################################################
            C = calc_measure(rnd_Mnet, old_linking_set, rnd_LS_df, "cLSsim", K_df=rnd_K_df)
            PR_C = calc_measure(rnd_Mnet, linking_set, rnd_LS_df, "new_cLS_PR", K_df=rnd_K_df)
            CinHubs_20 = calc_measure(rnd_Mnet, linking_set, rnd_LS_df, "CinLHubs_20", K_df=rnd_K_df)
            CinHubs_5 = calc_measure(rnd_Mnet, linking_set, rnd_LS_df, "CinLHubs_5", K_df=rnd_K_df)
            CinHubs = calc_measure(rnd_Mnet, linking_set, rnd_LS_df, "CinLHubs", K_df=rnd_K_df)
            PR_CH = calc_measure(rnd_Mnet, linking_set, rnd_LS_df, "PRofCinLhubs", K_df=rnd_K_df)
            #######################################################################
            df.loc[rep, ["name_set_a", "name_set_b", "linking_set", "name_layer_A", "name_layer_B"]] = [species_set_a, species_set_b, linking_set, name_layer_A, name_layer_B]
            df.loc[rep, ["N_A", "N_B", "K_A", "K_B", "C_A", "C_B", "N_lsA", "N_lsB"]] = \
            [N_A, N_B, K_A, K_B, C_A, C_B, N_lsA, N_lsB]
            df.loc[rep, ["Na", "Nb", "Nls"]] = [Na, Nb, Nls]
            df.loc[rep, ["HD", "LS_HD", "nonLS_HD", "P_HD", "nonP_HD", "lsA_HD", "lsB_HD", "r_k","rb_k"]] = [HD, LS_HD, nonLS_HD,P_HD,nonP_HD,ls_A_HD, ls_B_HD,r,rb]
            df.loc[rep, ["cLSsim","CinLHubs_20","CinLHubs_10","CinLHubs_5","new_cLS_PR","PRofCinLhubs"]]=[C,CinHubs_20,CinHubs,CinHubs_5,PR_C,PR_CH]

        df.to_csv(filename_df)
        #print(df.mean())
        #print(df.std())

        all_nets_av_data.loc[(sign, int, name),list(df)] = df.mean()
        all_nets_av_data.loc[(sign, int, name), ["name_set_a", "name_set_b", "linking_set", "name_layer_A", "name_layer_B"]] = [species_set_a,species_set_b,linking_set,name_layer_A,name_layer_B]
        all_nets_std_data.loc[(sign, int, name), list(df)] = df.std()
        all_nets_std_data.loc[(sign, int, name), ["name_set_a", "name_set_b", "linking_set", "name_layer_A", "name_layer_B"]] = [species_set_a, species_set_b, linking_set, name_layer_A, name_layer_B]

    #######################################################################
    # extinction AREA
    filename = "../OUTPUT/Data/Ext_Areas/Ext_Area_%s_%s_%s.csv" % (name, ext_MODE, null_model)
    rnd_Area_df = pd.read_csv(filename, index_col=[0, 1])
    if (sign == "AA"):  # correct column Area parasitism
        try:
            rnd_Area_df["Area_parasitism"] = rnd_Area_df["Area_Parasitoid"]
        except:
            rnd_Area_df["Area_parasitism"] = rnd_Area_df["Area_leaf_miner_parasitoid"]
    rnd_Area_df = rnd_Area_df.reorder_levels(["rnd", "ext"], axis=0)
    rnd_Area_df.sort_index(inplace=True)
    n_nets = rnd_Area_df.index.get_level_values("rnd").unique()
    n_exts = rnd_Area_df.index.get_level_values("ext").unique()
    print("%s nets ; %s ext sequences" % (len(n_nets),len(n_exts)))

    rnd_EA_mean = rnd_Area_df.mean()
    rnd_EA_std= rnd_Area_df.std()

    all_nets_av_data.loc[(sign, int, name), "Area_merged"] = rnd_EA_mean["Area_merged"]
    all_nets_av_data.loc[(sign, int, name), "EA_a"] = rnd_EA_mean.loc["Area_%s" % name_layer_A]
    all_nets_av_data.loc[(sign, int, name), "EA_b"] = rnd_EA_mean.loc["Area_%s" % name_layer_B]
    all_nets_av_data.loc[(sign, int, name), "EA_lA"] = rnd_EA_mean.loc["Area_only_%s" % name_layer_A]
    all_nets_av_data.loc[(sign, int, name), "EA_lB"] = rnd_EA_mean.loc["Area_only_%s" % name_layer_B]

    all_nets_std_data.loc[(sign, int, name), "Area_merged"] = rnd_EA_std["Area_merged"]
    all_nets_std_data.loc[(sign, int, name), "EA_a"] = rnd_EA_std.loc["Area_%s" % name_layer_A]
    all_nets_std_data.loc[(sign, int, name), "EA_b"] = rnd_EA_std.loc["Area_%s" % name_layer_B]
    all_nets_std_data.loc[(sign, int, name), "EA_lA"] = rnd_EA_std.loc["Area_only_%s" % name_layer_A]
    all_nets_std_data.loc[(sign, int, name), "EA_lB"] = rnd_EA_std.loc["Area_only_%s" % name_layer_B]

    all_nets_av_data.loc[(sign, int, name), "r_EA"] = rnd_Area_df.corr(method="spearman").loc["Area_%s" % name_layer_A, "Area_%s" % name_layer_B]
    corr_values=[]
    for n_net in n_nets:
        df = rnd_Area_df.xs(n_net, level="rnd")
        corr= df.corr(method="spearman").loc["Area_%s" % name_layer_A, "Area_%s" % name_layer_B]
        corr_values.append(corr)
    all_nets_std_data.loc[(sign, int, name), "r_EA"] = np.std(corr_values)
    #print("mean from each: %s " % np.mean(corr_values))
    #print("mean from all: %s" % rnd_Area_df.corr(method="spearman").loc["Area_%s" % name_layer_A, "Area_%s" % name_layer_B])
    #print("std:%s " % np.std(corr_values))

print("arrived at end, output information")
filename_df_av= "../OUTPUT/Data/Networks_rnd_df_%s_%s_%s_av.csv"  % (ext_MODE,null_model,Nrep)
print(filename_df_av)
all_nets_av_data.to_csv(filename_df_av)
filename_df_std= "../OUTPUT/Data/Networks_rnd_df_%s_%s_%s_std.csv" % (ext_MODE,null_model,Nrep)
all_nets_std_data.to_csv(filename_df_std)
print(filename_df_std)
