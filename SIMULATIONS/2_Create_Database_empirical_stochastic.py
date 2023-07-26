from Mnetworks import *
from Metrics import *
import sys


#1st decide the networks you want to study: in this case only MA 
Mnet_sign_types=["MA"] #"MMM" , "MAA"
# #
network_names=[]
for net_sign in Mnet_sign_types:
     #print("net sign:%s" % net_sign)
     for net_type in dict_sign_types_2_net_type[net_sign]:
        #if(net_type == "MA" ):
        #print("net type:%s\n" % net_type)
        network_names.append(dict_type_2_names[net_type])
network_names=network_names[0]

network_name_type_tuples=[]
for name in network_names:
     network_name_type_tuples.append((dict_net_type[dict_name_net[name]], dict_name_net[name], name))
#create index
my_index = pd.MultiIndex.from_tuples(network_name_type_tuples,names=['sign','int','name'])

#define the extinction scenario we want to study
ext_MODE="RND" #, DD, ID

pd.options.mode.chained_assignment = None  # default='warn'

#determine the structural metrics we want to include in the database
measures=["name_set_a","name_set_b","linking_set","name_layer_A","name_layer_B","Na","Nb","Nls","N_lsA","N_lsB","N_A","N_B","K_A","K_B","C_A","C_B", "HD","LS_HD","nonLS_HD", "P_HD", "nonP_HD","lsA_HD","lsB_HD" ,"rb_k", "r_k","cLSsim","CinLHubs_10","CinLHubs_20","CinLHubs_5","new_cLS_PR","PRofCinLhubs", "Area_merged", "Area_merged_std","EA_a","EA_b","EA_lA","EA_lB","r_EA","AreaW_merged","AreaW_merged_std","EAW_a","EAW_b","EAW_lA","EAW_lB","r_EAW"]


data=pd.DataFrame(index=my_index, columns=measures)

for sign, int, name in network_name_type_tuples:
    print(name)
    Mnet=Read_net_general(name)
    linking_set=get_linking_set(Mnet)
    K_df, linking_set = get_full_K_df(name,Mnet,linking_set=linking_set)
    LS_df=K_df[K_df["set"]==linking_set]
    old_linking_set = get_linking_set(Mnet)
    #size and connectance related stuff #########################################################
    #get size of subnetworks
    sub_net_sizes={}
    sub_net_links={}
    max_size=0
    max_name=0
    interactions = list(Mnet.slices[2])
    for interaction in interactions:
        nodes, links =get_nodes_in_layer(Mnet,interaction_layer=interaction, nlinks=True)
        size=len(nodes)
        sub_net_sizes[interaction]=size
        sub_net_links[interaction]=links
        if (size>max_size):
            max_size=size
            max_name=interaction

    name_layer_A = max_name
    N_A=sub_net_sizes.pop(name_layer_A)
    L_A = sub_net_links.pop(name_layer_A)
    K_A=L_A/N_A

    name_layer_B=list(sub_net_sizes.keys())[0]
    N_B=sub_net_sizes.pop(name_layer_B)
    L_B = sub_net_links.pop(name_layer_B)
    K_B=L_B/N_B


    #get sizes of sets
    sp_sets=list(K_df["set"].unique())
    #get sizes os linking set species in different layers
    nodes_lsA=get_nodes_in_spset(Mnet,sign,species_set= linking_set,interaction_layer=name_layer_A)
    N_lsA=len(nodes_lsA)
    nodes_lsB=get_nodes_in_spset(Mnet,sign,species_set=linking_set,interaction_layer=name_layer_B)
    N_lsB=len(nodes_lsB)
    Nls=len(set(get_nodes_in_spset(Mnet,sign,species_set=linking_set)))

    species_set_a=dict_sp_set_name_from_interaction[sign][name_layer_A]
    species_set_b=dict_sp_set_name_from_interaction[sign][name_layer_B]
    Na=len(get_nodes_in_spset(Mnet, sign,species_set=species_set_a ,interaction_layer=name_layer_A))
    Nb=len(get_nodes_in_spset(Mnet, sign,species_set=species_set_b  , interaction_layer=name_layer_B))

    #get connectances of each layer:
    C_A=(L_A*0.5)/(Na*N_lsA)
    C_B=(L_B*0.5)/(Nb*N_lsB)

    # fill basic structural measure
    # "EA", "EA_a","EA_b","EA_lA","EA_lB","r_EA"
    # by layer
    data.loc[(sign, int, name), ["name_set_a", "name_set_b", "linking_set", "name_layer_A", "name_layer_B"]] = [species_set_a, species_set_b, linking_set, name_layer_A, name_layer_B]
    data.loc[(sign, int, name), ["N_A", "N_B", "K_A", "K_B", "C_A", "C_B", "N_lsA", "N_lsB"]] = [N_A, N_B, K_A, K_B, C_A, C_B, N_lsA, N_lsB]
    #by set
    data.loc[(sign, int, name),["Na","Nb","Nls"]]=[Na, Nb, Nls]

    # Degree distribution related stuff ######################################################
    #Degree heterogeneities
    HD=calc_DH(K_df)
    LS_HD=calc_DH(LS_df)
    nonLS_HD=calc_DH(K_df[K_df["set"]!=linking_set])

    P_HD=calc_DH(K_df[K_df["set"]=="Plant"])
    nonP_HD=calc_DH(K_df[K_df["set"]!="Plant"])

    ls_A_HD = (K_df.loc[nodes_lsA].std() / K_df.loc[nodes_lsA].mean())[name_layer_A]
    ls_B_HD = (K_df.loc[nodes_lsB].std() / K_df.loc[nodes_lsB].mean())[name_layer_B]

    #Degree-degree correlations
    rb = calc_measure(Mnet,old_linking_set, LS_df, "rb_k",K_df=K_df)
    r = calc_measure(Mnet,old_linking_set,LS_df,"r_k",K_df=K_df)


    #fill degree heterogeneities and degree-degree correlations
    #print([HD, LS_HD,nonLS_HD,P_HD,nonP_HD,ls_A_HD,ls_B_HD,rb])
    #print(data.loc[(sign, int, name),["HD","LS_HD","nonLS_HD", "P_HD", "nonP_HD","lsA_HD","lsB_HD", "rb_k","r_k"]])
    data.loc[(sign, int, name),["HD","LS_HD","nonLS_HD", "P_HD", "nonP_HD","lsA_HD","lsB_HD", "rb_k","r_k"]]=[HD, LS_HD,nonLS_HD,P_HD,nonP_HD,ls_A_HD,ls_B_HD,rb,r]

    # Linking set metrics ###################################################################
    C=calc_measure(Mnet,old_linking_set,LS_df,"cLSsim",K_df=K_df)
    PR_C = calc_measure(Mnet, linking_set, LS_df, "new_cLS_PR", K_df=K_df)
    CinHubs_20=calc_measure(Mnet,linking_set,LS_df,"CinLHubs_20",K_df=K_df)
    CinHubs_5 = calc_measure(Mnet, linking_set, LS_df, "CinLHubs_5", K_df=K_df)
    CinHubs=calc_measure(Mnet,linking_set,LS_df,"CinLHubs",K_df=K_df)
    PR_CH=calc_measure(Mnet,linking_set,LS_df,"PRofCinLhubs",K_df=K_df)


    data.loc[(sign, int, name),["cLSsim","CinLHubs_20","CinLHubs_10","CinLHubs_5","new_cLS_PR","PRofCinLhubs"]]=[C,CinHubs_20,CinHubs,CinHubs_5,PR_C,PR_CH]

    # Extinction areas
    # #"EA", "EA_a","EA_b","EA_lA","EA_lB","r_EA"]#####################################################################
    #EA_filename= "./OUTPUT/Data/Ext_Area_%s_%s.csv" % (name, ext_MODE)
    EA_filename= "../../OUTPUT/Data/Ext_Area_%s_%s.csv" % (name, ext_MODE)
    
    EA_df=pd.read_csv(EA_filename,index_col=0)
    if (sign == "AA"):#correct column Area parasitism
        try:
            EA_df["Area_parasitism"]=EA_df["Area_Parasitoid"]
        except:
            EA_df["Area_parasitism"]=EA_df["Area_leaf_miner_parasitoid"]
    mean_EA=EA_df.mean()
    std_EA = EA_df.std()

    data.loc[(sign, int, name),"Area_merged"]=mean_EA["Area_merged"]
    data.loc[(sign, int, name), "Area_merged_std"] = std_EA["Area_merged"]
    data.loc[(sign, int, name),"EA_a"]=mean_EA["Area_%s" % name_layer_A]
    data.loc[(sign, int, name), "EA_b"] = mean_EA["Area_%s" % name_layer_B]
    data.loc[(sign, int, name), "EA_lA"] = mean_EA["Area_only_%s" % name_layer_A]
    data.loc[(sign, int, name), "EA_lB"] = mean_EA["Area_only_%s" % name_layer_B]

    data.loc[(sign, int, name), "r_EA"] = EA_df.corr(method="spearman").loc["Area_%s" % name_layer_A, "Area_%s" % name_layer_B]

    #weighted extinction areas
    #EAW_filename= "./OUTPUT/Data/Ext_Area_%s_%s_W.csv" % (name, ext_MODE)
    EAW_filename= "../../OUTPUT/Data/Ext_Area_%s_%s_W.csv" % (name, ext_MODE)
    EAW_df=pd.read_csv(EAW_filename,index_col=0)
    mean_EAW=EAW_df.mean()
    std_EAW = EAW_df.std()

    data.loc[(sign, int, name),"AreaW_merged"]=mean_EAW["Area_merged"]
    data.loc[(sign, int, name), "AreaW_merged_std"] = std_EAW["Area_merged"]
    data.loc[(sign, int, name),"EAW_a"]=mean_EAW["Area_%s" % name_layer_A]
    data.loc[(sign, int, name), "EAW_b"] = mean_EAW["Area_%s" % name_layer_B]
    data.loc[(sign, int, name), "EAW_lA"] = mean_EAW["Area_only_%s" % name_layer_A]
    data.loc[(sign, int, name), "EAW_lB"] = mean_EAW["Area_only_%s" % name_layer_B]

    data.loc[(sign, int, name), "r_EAW"] = EAW_df.corr(method="spearman").loc["Area_%s" % name_layer_A, "Area_%s" % name_layer_B]


#print(data)
df_filename="../OUTPUT/Data/Networks_df_%s_W.csv" % ext_MODE
#print(df_filename)
data.reset_index(inplace=True)
data.to_csv(df_filename)

df=pd.read_csv(df_filename, index_col=0)