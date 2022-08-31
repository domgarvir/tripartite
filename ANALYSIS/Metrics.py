import numpy as np

from Mnetworks import *

# CALCULATING PROPERTIES
#measure Participation Ratio of nodes
def calc_new_PR(K_df,get_linking_set):
    cols_to_use = list(set(list(K_df)).intersection(set(positive_interactions).union(negative_interactions)))

    #now PR is just the minimum relative degree *2
    modDfObj = K_df[cols_to_use].div(K_df["O.Degree"], axis=0)
    new_PR=modDfObj.min(axis=1)*2
    return new_PR

#measure average degree of neighbours of nodes
def calc_Knn(K_df,Mnet,**kwargs):
    interactions=list(Mnet.slices[2])
    for nodei in list(Mnet._net.keys()):
        ki=len(Mnet._net[nodei].keys())
        #print("%s,%s,%s (%s) se conecta con:" % (nodei[0],nodei[1],nodei[2],ki))
        knn={interactions[0]:0, interactions[1]:0, "merged":0}
        k={interactions[0]:0, interactions[1]:0, "merged":0}
        # for nodej in list(Mnet._net[nodei].keys()):
        #     kj=len(Mnet._net[nodej].keys())
        #     print("%s,%s,%s (%s)" % (nodej[0],nodej[1],nodej[2],kj))
        #     knn[nodei[2]] = knn[nodei[2]] + kj
        for nodej in list(Mnet._net[nodei].keys()):

            for interaction in interactions:
                try:
                    kj=len(Mnet._net[(nodej[0],nodej[1],interaction)].keys())
                    #print("%s,%s,%s (%s)" % (nodej[0], nodej[1], interaction, kj))
                    knn[interaction] = knn[interaction] + kj
                    k[interaction] = k[interaction] + 1
                    knn["merged"] = knn["merged"] + kj
                    k["merged"]= k["merged"] +1
                except:
                    pass
        for interaction in interactions:
            if (k[interaction]>0):
                try:
                    K_df.loc[nodei[0], "Knn_%s" % interaction] = knn[interaction] / k[interaction]
                    ##K_df.loc[nodei[0], "Knn_%s" % interaction] = knn[interaction]
                    K_df.loc[nodei[0], "K_%s" % interaction] = k[interaction]
                except:
                    pass

            else:
                try:
                    if (K_df.loc[nodei[0],"O.Degree"] == K_df.loc[nodei[0],"%s" % interactions[0]]):
                        K_df.loc[nodei[0], "Knn_%s" % interaction] = 0
                    #else:
                except:
                    pass

        if(nodei[0] in list(K_df.index)):
            K_df.loc[nodei[0],"Knn_merged"] = knn["merged"]/k["merged"]

    K_df.fillna(0,inplace=True)

    try:
        only_Knn = kwargs["only_Knn"]
        if (only_Knn == True):
            for interaction in interactions:
                K_df.drop("K_%s" % interaction, axis=1, inplace=True)

    except:
        pass

    return K_df

#measure average Participation ratio of node neighbours
def calc_PRnn(K_df,Mnet,**kwargs):

    linking_set=get_linking_set(Mnet)
    interactions = list(Mnet.slices[2])

    PR_metric="new_PR"
    try:
        PR_metric=kwargs["PR_metric"]
    except:
        pass

    K_df["%s_nn_%s" % (PR_metric,interactions[0])]=0
    K_df["%s_nn_%s" % (PR_metric,interactions[1])] = 0

    PR_nn = {interactions[0]: 0, interactions[1]: 0, "merged": 0}
    non_LS_nodes=K_df[K_df["set"]!=linking_set].index.tolist()

    for nodei in list(Mnet._net.keys()):#for each non linking set node, look at the average PR of their neighbours
      PR_nn = {interactions[0]: 0, interactions[1]: 0, "merged": 0}
      if (nodei[1] != linking_set):
        #print("\nnodei:%s in %s" % (nodei[0],nodei[2]))
        interaction = nodei[2]
        for nodej in list(Mnet._net[nodei].keys()):#for each partner
            try:
                PR_nn_i=K_df.loc[nodej[0],PR_metric]
                PR_nn[interaction] =PR_nn[interaction]+PR_nn_i
                #print("%s in %s, PR:%s" % (nodej[0], nodej[2], PR_nn_i))
            except:
                #print("%s not found!" %nodej[0])
                pass
        #print("PR_nn total =%s , divided by %s gives %s\n" % (PR_nn[interaction],max(1,len(Mnet._net[nodei].keys())),PR_nn[interaction]/(max(1,len(Mnet._net[nodei].keys())))))
        PR_nn[interaction]=PR_nn[interaction]/(max(1,len(Mnet._net[nodei].keys())))
        K_df.loc[nodei[0],"%s_nn_%s" % (PR_metric,interaction)]=PR_nn[interaction]

    K_df["%s_nn_merged" % PR_metric]=K_df["%s_nn_%s" % (PR_metric,interactions[0])] + K_df["%s_nn_%s" % (PR_metric,interactions[1])]

    return K_df

#measure degree heterogeneity
def calc_DH(K_df):
    av_sigma_LS=(K_df.std()/K_df.mean())["O.Degree"]
    return av_sigma_LS

def calc_CinHubs(K_df,treshold=10,measure="CinH"):

    df=K_df.sort_values(by=["O.Degree"], ascending=False)

    #subselect 10% (or X%) more connected
    lastN=max(1,int(K_df.shape[0]*treshold/100))
    #get degree of last node
    last_k = df.iloc[lastN-1]["O.Degree"]
    # now the average PR for all those with that degree
    av_PR = df.groupby("O.Degree").mean().loc[last_k, "new_PR"]
    #nb of nodes with same k as last
    nb_of_nodes_k=df.groupby("O.Degree").count().loc[last_k,"new_PR"]
    prob_of_connectors_with_k=((df[df["O.Degree"]==last_k]["new_PR"]!=0).sum())/nb_of_nodes_k

    df_hubs=df.iloc[0:lastN,:]
    nb_of_nodes_k_in_hubs= df_hubs.groupby("O.Degree").count().loc[last_k, "new_PR"]

    if(measure=="CinH"):
        CinH=(df_hubs[df_hubs["O.Degree"]>last_k]['new_PR'] != 0).sum()
        CinH += nb_of_nodes_k_in_hubs*prob_of_connectors_with_k
        return CinH / lastN
    elif(measure=="PRofCinH"):
        PRofCinH=df_hubs[df_hubs['new_PR']>0].mean()["new_PR"]
        if (np.isnan(PRofCinH)):
            PRofCinH=0
        return PRofCinH
    elif (measure == "PRinH"):
        PRofCinH = df_hubs['new_PR'].mean()
        return PRofCinH
    elif(measure=="IPRofCinLhubs"):
        IPRofCinLhubs=df_hubs['IPR'].mean()
        return IPRofCinLhubs
    elif (measure=="ALL"):
        CinH = (df_hubs['new_PR'] != 0).sum()
        PRofCinH = df_hubs['new_PR'].mean()
        return CinH / lastN, PRofCinH

def calc_degree_degree_correlations_bipart(Mnet,linking_set):

    #we will have two sets of nodes: linking nodes [A] and the others [B]
    A=get_full_nodes_in_set(Mnet,set_name=linking_set)
    B=[]
    for eachset in list(Mnet.slices[1]):
        if(eachset != linking_set):
            B.append(get_full_nodes_in_set(Mnet,set_name=eachset))
    B=list(itertools.chain.from_iterable(B))
    N={"A":0,"B":0}
    k_med={"A":0,"B":0}
    k2_med={"A":0,"B":0}
    k3_med={"A":0,"B":0}

    node_sets={"A":A,"B":B}

    for each_set in node_sets.keys():
        for node in node_sets[each_set]:
            k=len(Mnet._net[node].keys())
            k_med[each_set] = k_med[each_set] + k
            k2_med[each_set] = k2_med[each_set] + pow(k, 2)
            k3_med[each_set] = k3_med[each_set] + pow(k, 3)

        N[each_set] = len(node_sets[each_set])
        k_med[each_set]  =  k_med[each_set] / N[each_set]
        k2_med[each_set] = k2_med[each_set] / N[each_set]
        k3_med[each_set] = k3_med[each_set] / N[each_set]

    # print("-------")
    # print(k_med,k2_med,k3_med,N,node_sets)
    # print("-------")

    L=0
    kikj_med=0
    for nodei in A:
        ki=len(Mnet._net[nodei].keys())
        for nodej in Mnet._net[nodei].keys():
            kj=len(Mnet._net[nodej].keys())
            kikj_med = kikj_med + ki*kj

            L = L + 1

    kikj_med = kikj_med / L
    #print(kikj_med)

    kmed2=(k2_med["A"]*N["A"])/(k_med["A"]*N["B"])
    kmed1=(k2_med["B"]*N["B"])/(k_med["B"]*N["A"])
    k2med2=(k3_med["A"]*N["A"])/(k_med["A"]*N["B"])
    k2med1=(k3_med["B"]*N["B"])/(k_med["B"]*N["A"])


    #r_pearson= ( kikj_med - (pow(k2_med,2)/pow(k_med,2)) )/( (k3_med/k_med) - (pow(k2_med,2)/pow(k_med,2)))

    sigma1 = (((k3_med["A"] * N["A"]) / (k_med["A"] * N["B"])) - (((k2_med["A"] * N["A"]) / (k_med["A"] * N["B"])) * ((k2_med["A"] * N["A"]) / (k_med["A"] * N["B"])))); # Teoria 1

    sigma2 = (((k3_med["B"] * N["B"]) / (k_med["B"] * N["A"])) - (((k2_med["B"] * N["B"]) / (k_med["B"] * N["A"])) * ((k2_med["B"] * N["B"]) / (k_med["B"] * N["A"]))));

    #print(kmed2,kmed1,sigma1,sigma2)

    r_pearson = (kikj_med - (kmed1 * kmed2)) / (np.sqrt(np.fabs(sigma1 * sigma2)));


    return r_pearson, k_med, k2_med, k3_med
def calc_degree_degree_correlations(Mnet):

    r_pearson=0
    k_med=0
    k2_med=0
    k3_med=0
    N=len(Mnet._net.keys())
    for node in Mnet._net.keys():
        k=len(Mnet._net[node].keys())
        k_med = k_med + k
        k2_med = k2_med + pow(k,2)
        k3_med = k3_med + pow(k,3)

    k_med = k_med/N
    k2_med = k2_med/N
    k3_med = k3_med/N

    L=0
    kikj_med=0
    for nodei in Mnet._net.keys():
        ki=len(Mnet._net[nodei].keys())
        for nodej in Mnet._net[nodei].keys():
            kj=len(Mnet._net[nodej].keys())
            kikj_med = kikj_med + ki*kj

            L = L + 1

    kikj_med = kikj_med / L
    #print(kikj_med)
    r_pearson= ( kikj_med - (pow(k2_med,2)/pow(k_med,2)) )/( (k3_med/k_med) - (pow(k2_med,2)/pow(k_med,2)))

    return r_pearson, k_med, k2_med, k3_med

def calc_measure(Mnet,linking_set,LS_df,measure,**kwargs):
    K_df = None
    try:
        K_df = kwargs["K_df"]
    except:
        pass

    if (measure == "cLSsim"):
        connectors=LS_df[LS_df["new_PR"]!=0].count()[0]
        all=LS_df.shape[0]
        return connectors/all

    elif (measure=="new_cLS_PR") :#"measures of PR in connector nodes")
        LS_df_clean=LS_df[LS_df["new_PR"]!=0]
        return LS_df_clean["new_PR"].mean()

    elif (measure == "new_PR"): #measures of PR in linking set nodes
        return LS_df["new_PR"].mean()

    elif ((measure.startswith("CinLHubs"))): #measure prop of Connectors in Hubs
        try:
            m,treshold=measure.rstrip().split("_")
        except:
            treshold=10
        CinLHubs=calc_CinHubs(LS_df,treshold=int(treshold))
        return CinLHubs

    elif (measure=="PRofCinLhubs"):# Participation rato of the hub-connectors
        PRofCinLhubs=calc_CinHubs(LS_df,treshold=20,measure="PRofCinH")
        return PRofCinLhubs

    elif (measure == "rb_k"):#degree-degree correlation in bipartite network
        MMnet = get_merged_multinetwork(Mnet)
        r_all = calc_degree_degree_correlations_bipart(MMnet,linking_set)
        return r_all[0]

    elif (measure == "r_k"):
        MMnet = get_merged_multinetwork(Mnet)
        r_all = calc_degree_degree_correlations(MMnet)
        return r_all[0]
    else:
        print("ERROR %s measure not found, return None" % measure)

# get basic dataframe of connectivity and nodes properties
def get_K_df(Mnet,linking_set,marker="_",**kwargs):

    aspects = list(Mnet.slices[2]) #different networks
    sets = list(Mnet.slices[1]) #different groups of nodes (plants, pollinators...)

    #cretae all the column names for the dataframe
    cols=["set"]
    for aspect in aspects:
        cols.append(aspect)


    #create index for the dataframe and initialize with 0s
    nodes=[]
    for key in Mnet._net.keys():
        nodes.append(key[0])
    nodes=list(set(nodes))
    #initialize dataframe of konnectivity. Species in rows, different networks in columns, incliding species type
    K_df=pd.DataFrame(0,index=nodes,columns=cols)
    double_species_df=pd.DataFrame(columns=cols)
    sp_to_merge=[]
    #print(K_df)
    for key in Mnet._net.keys():
        if (K_df.loc[key[0],"set"]==0):
            K_df.loc[key[0],"set"]=key[1]
            K_df.loc[key[0], key[2]] = len(list(Mnet._net[key]))
        else:
            if(K_df.loc[key[0], "set"] != linking_set ):
                if (K_df.loc[key[0], key[2]]==0):
                    #K_df.loc[key[0], "set"] = K_df.loc[key[0], "set"] + "&" + key[1]
                    di={'name':"%s%s%s" % (key[0],marker,key[2]),'set': key[1], key[2]: len(list(Mnet._net[key]))}
                    double_species_df=double_species_df.append(di,ignore_index=True)
                else:
                    K_df.loc[key[0], key[2]] = K_df.loc[key[0], key[2]]+ len(list(Mnet._net[key]))

                    if (key[1] != linking_set):
                        #print("añadimos esta como doble:")
                        #print(key)
                        sp_to_merge.append((key[0],key[2]))
        #print("%s: set %s %s" % (key[0],key[1],K_df.loc[key[0],"set"]))
            else:
                K_df.loc[key[0], key[2]] = len(list(Mnet._net[key]))

    #now for each species in two sets, duplicate it

    try: #try to add the double set species (if there are)
        double_species_df.set_index('name',inplace=True)
        double_species_df.fillna(0, inplace=True)
        K_df=K_df.append(double_species_df)
    except:
        pass



    cols_to_use=list(set(list(cols)).intersection(set(positive_interactions).union(negative_interactions)))
    K_df.loc[:, "O.Degree"] = K_df[cols_to_use].sum(axis=1)
    K_df.loc[:,"new_PR"]=calc_new_PR(K_df,linking_set)
    K_df.fillna(0,inplace=True)

    #print(K_df)
    return K_df, double_species_df, sp_to_merge

#get full dataframe of nodes properties
def get_full_K_df(name,Mnet, marker="_",**kwargs):

    return_double_set_species=False
    try:
        return_double_set_species=kwargs["return_2set_sp"]
    except:
        pass


    linking_set="Plant"
    try:
        linking_set=kwargs["linking_set"]
    except:
        pass

    K_df , K_2set_df, sp_to_merge = get_K_df(Mnet, marker=marker,linking_set=linking_set)
    int_type = dict_name_net[name]
    sign_type = dict_net_type[int_type]
    K_df.replace({"set": dict_node_simple[sign_type]},inplace=True)
    linking_set = dict_node_simple[sign_type][linking_set]
    K_df.index.name = "sp_name"
    K_df = calc_Knn(K_df, Mnet)
    K_df = calc_PRnn(K_df, Mnet, PR_metric="new_PR")
    K_df.loc[:,"O.D_zsc"] = (K_df.loc[:,"O.Degree"] - K_df.loc[:,"O.Degree"].mean()) / (K_df.loc[:,"O.Degree"].std())
    # print(K_df[K_df["O.D_zsc"]>1.65])
    K_df.loc[:,"Knn_zsc"] = (K_df['Knn_merged'] - K_df['Knn_merged'].mean()) / (K_df['Knn_merged'].std())

    if (return_double_set_species):
        return K_df, linking_set, K_2set_df, sp_to_merge
    else:
        return K_df, linking_set

def calc_Z(ma, mb, sa, sb,nb_samples=3000):
    Z = (ma - mb) / np.sqrt(pow(sa/np.sqrt(nb_samples), 2) + pow(sb/np.sqrt(nb_samples), 2))
    return Z

def get_ranking_df():
    df_rank_similitude = pd.DataFrame()
    index = 0

    filtered_network_names = ['Sinohara_1_ALL_PH', 'Sinohara_2_ALL_PH', 'Sinohara_3_ALL_PH', 'Sinohara_4_ALL_PH', 'Sinohara_ALL_A_PH', 'Sinohara_ALL_E_PH', 'Sinohara_ALL_I_PH', 'Sinohara_2_E_PH','Sinohara_3_E_PH', 'Sinohara_4_I_PH', 'Melian_OO_OO_PH', 'Hackett_1_ALL_PH', 'Hackett_2_ALL_PH', 'Hackett_1_S_PH', 'Hackett_1_GL_PH', 'Pocock_OO_OO_PH','Melian_OO_OO_HSD', 'McFayden_ALL_A_HPa', 'McFayden_1_A_HPa', 'McFayden_2_A_HPa','McFayden_3_A_HPa', 'McFayden_4_A_HPa', 'McFayden_5_A_HPa', 'McFayden_6_A_HPa', 'McFayden_7_A_HPa', 'McFayden_8_A_HPa', 'McFayden_9_A_HPa', 'McFayden_10_A_HPa', 'McFayden_ALL_B_HPa', 'McFayden_1_B_HPa', 'McFayden_2_B_HPa', 'McFayden_3_B_HPa',  'McFayden_4_B_HPa', 'McFayden_5_B_HPa', 'McFayden_6_B_HPa', 'McFayden_7_B_HPa','McFayden_8_B_HPa', 'McFayden_9_B_HPa', 'McFayden_10_B_HPa', 'Hackett_1_ALL_HPa','Hackett_1_WL_HPa', 'Melian_OO_OO_PSD', 'Dattilo_OO_OO_PSD', 'Dattilo_OO_OO_PA']

    df_filename = "../OUTPUT/Data/Networks_df_%s.csv" % "RND"
    df = pd.read_csv(df_filename, header=0, index_col="name")

    for name in filtered_network_names:
        print(name)
        my_sign = dict_net_type[dict_name_net[name]]
        print(my_sign)
        Mnet = Read_net_general(name)
        linking_set = get_linking_set(Mnet)
        nodes_to_erase = get_nodes_in_set(Mnet, set_name="Plant")
        N = len(nodes_to_erase)
        interactions = list(Mnet.slices[2])
        K_df, linking_set, trash = get_K_df(Mnet, linking_set=linking_set)

        filename =  "../OUTPUT/Data/Ext_Areas/Node_impact_%s_%s.csv" % (name, "RND")
        Node_impact_df = pd.read_csv(filename, index_col=0)
        # print(Node_impact_df.mean()[["Area_%s" % interactions[0],"Area_%s" % interactions[1]]])
        if (my_sign == "AA"):  # correct the error in Area parsitism
            # print("CHANGE!###############################################")
            try: #unifying names of interaction
                Node_impact_df["Area_parasitism"] = Node_impact_df["Area_Parasitoid"]
            except:
                Node_impact_df["Area_parasitism"] = Node_impact_df["Area_leaf_miner_parasitoid"]
        # print(Node_impact_df.mean()[["Area_%s" % interactions[0], "Area_%s" % interactions[1]]])

        Node_corr_df = Node_impact_df.corr(method="spearman")

        rankings_names = ["Area_merged", "Area_%s" % interactions[0], "Area_%s" % interactions[1],
                          "Area_only_%s" % interactions[0], "Area_only_%s" % interactions[1]]
        ranking = {}

        for rname in rankings_names:
            ranking[rname] = Node_corr_df[rname][:N].rank(ascending=False)
            ranking[rname].name = rname
            ranking[rname].to_frame().reset_index()

        list_of_rankings = list(ranking.values())
        all_ranks = pd.concat(list_of_rankings, axis=1)
        all_ranks = all_ranks.sort_values("Area_merged", ascending=False)
        #print(all_ranks.corr(method="spearman"))


        columns = ["Area_%s" % interactions[0], "Area_merged", "Area_%s" % interactions[1]]

        r_corr_df = all_ranks[columns].corr(method="spearman")
        keep = np.triu(np.ones(r_corr_df.shape)).astype('bool').reshape(r_corr_df.size)
        r_corr_vertical = r_corr_df.stack()[keep].reset_index()
        r_corr_vertical["r_EA"] = df.loc[name]["r_EA"]
        r_corr_vertical["name"] = name
        r_corr_vertical["sign"] = dict_net_type[dict_name_net[name]]

        df_rank_similitude = pd.concat([df_rank_similitude, r_corr_vertical])
        #print(r_corr_vertical)

    df_rank_similitude.rename(columns={0:'corr'},inplace=True)
    filename = "../OUTPUT/Data/All_Ranking_%s.csv" % ("RND")
    df_rank_similitude.to_csv(filename)
    print(df_rank_similitude)
    return df_rank_similitude
