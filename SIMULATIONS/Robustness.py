from Mnetworks import *
import copy
#create or initialize structures
def create_area_struct(Mnet,active_set,nodes_to_erase):
    nxMnet = from_pymnet_2_nx(Mnet)
    Area_struct = {}
    Area_struct["N_all"] = 0.
    for my_set in Mnet.slices[1]:
        Area_struct["N_%s" % my_set] = len(get_nodes_in_set(Mnet, set_name=my_set))
        Area_struct["N_all"] = Area_struct["N_all"] + Area_struct["N_%s" % my_set]
        Area_struct["Area_%s" % my_set] = 0.

    Area_struct["N_passive"] = Area_struct["N_all"] - Area_struct["N_%s" % active_set]
    Area_struct["N_active"] = Area_struct["N_%s" % active_set]
    Area_struct["Area_merged"] = 0.

    for interaction in Mnet.slices[2]:
        nodes_in_net = set(nxMnet[interaction].nodes)
        Area_struct["N_active_%s" % interaction] = len(nodes_in_net & set(nodes_to_erase))
        Area_struct["N_passive_%s" % interaction] = len(nodes_in_net - set(nodes_to_erase))
        Area_struct["Area_%s" % interaction] = 0.
        Area_struct["Area_only_%s" % interaction] = 0.

    return Area_struct
def create_ext_struct(Mnet):
    ext_struct = {}
    ext_struct["merged"] = 0
    for my_set in Mnet.slices[1]:
        ext_struct[my_set] = 0

    for interaction in Mnet.slices[2]:
        ext_struct[interaction] = 0
        ext_struct["Plant_%s" % interaction] = 0

    return ext_struct
def Initialize_Area_df(Mnet,**kwargs):

    sets=Mnet.slices[1]
    interactions=Mnet.slices[2]

    col_names=["Area_merged"]

    for each_set in sets:
        col_names.append("Area_%s" % each_set)

    for interaction in interactions:
        col_names.append("Area_%s" % interaction)
        col_names.append("Area_only_%s" % interaction)

    try :
        Nseq=kwargs["Nrep"]
    except:
        Nseq=False

    try:
        Nnet=kwargs["Nnet"]
    except:
        Nnet=False

    if ((Nnet!=0) & (Nseq!=0)):
        arrays=[np.arange(1,Nseq),np.arange(1,Nnet)]
        mi=pd.MultiIndex.from_product([arrays[1], arrays[0]],names=['rnd', 'ext'])
        area_df = pd.DataFrame(0,index=mi, columns=col_names)
    elif (Nseq !=0):
        area_df = pd.DataFrame(0,index=np.arange(Nseq), columns=col_names)
    else:
        area_df = pd.DataFrame(0,columns=col_names)
    # try:
    #     area_df=pd.DataFrame(index= np.arange(kwargs["Nrep"]),columns=col_names)
    # except:
    #     area_df = pd.DataFrame(columns=col_names)

    return area_df
def initialize_node_impact_df(Mnet,active_set,nodes_to_erase,area_struct,**kwargs):

    col_names=nodes_to_erase.copy()
    #col_names.append("Area_merged")

    sets = Mnet.slices[1]
    interactions = Mnet.slices[2]
    # for each_set in sets:
    #     col_names.append("Area_%s" % each_set)
    #
    # for interaction in interactions:
    #     col_names.append("Area_%s" % interaction)
    #     col_names.append("Area_only_%s" % interaction)

    for keys in area_struct.keys():
        col_names.append(keys)

    try:
        node_df=pd.DataFrame(index= np.arange(kwargs["Nrep"]),columns=col_names)
    except:
        node_df = pd.DataFrame(columns=col_names)

    return node_df
def reset_area_struct(area_struct,Mnet):

    area_struct["Area_merged"]=0.

    for my_set in Mnet.slices[1]:
        area_struct["Area_%s" % my_set]=0.

    for interaction in Mnet.slices[2]:
        area_struct["Area_%s" % interaction]=0.
        area_struct["Area_only_%s" % interaction] = 0.

    return area_struct
def reset_ext_struct(ext_struct,Mnet):

    ext_struct["merged"]=0
    for my_set in Mnet.slices[1]:
        ext_struct[my_set] = 0
    for interaction in Mnet.slices[2]:
        ext_struct[interaction] = 0
        ext_struct["Plant_%s" % interaction] = 0

    return ext_struct
# Add data to estructures
def Add_ext_area_to_Area_df(area_df,area_struct,seq,**kwargs):

    try:
        net=kwargs["net"]
    except:
        net=False


    ext_areas=pd.Series(0,index=list(area_df))

    for i in ext_areas.index:
        ext_areas.loc[i]=area_struct[i]

    if net:
        area_df.loc[(net,seq)] = ext_areas
    else:
        area_df.iloc[seq] = ext_areas
    return area_df

def Add_node_seq_to_node_impact_df(Node_impact_df,node_sequence,ext_area,rep):

    index=0
    for index in range(len(node_sequence)):
        Node_impact_df.loc[rep,node_sequence[index]]=index

    for akey in ext_area.keys() :
        Node_impact_df.loc[rep,akey]=ext_area[akey]


    return Node_impact_df
#########################################
def get_node_sequence(nodes_to_erase, linking_set_DF, ext_MODE):

    if (ext_MODE == "RND"):
        node_sequence = sample(nodes_to_erase, len(nodes_to_erase))
    else:
        cols_to_use=list(set(list(linking_set_DF)).intersection(set(positive_interactions).union(negative_interactions)))
        linking_set_DF.loc[:, "O.Degree"] = linking_set_DF[cols_to_use].sum(axis=1)
        node_sequence = []
        if(ext_MODE == "DD"):#decreasing degree
            raw_ranking = linking_set_DF.sort_values(by="O.Degree", ascending=False)["O.Degree"]
        elif (ext_MODE =="ID") :#increasing degree
            raw_ranking = linking_set_DF.sort_values(by="O.Degree", ascending=True)["O.Degree"]

        ranking_values = raw_ranking.unique()
        for value in ranking_values:
            tie_nodes = list(raw_ranking[raw_ranking == value].index)
            rnd_untie = sample(tie_nodes, len(tie_nodes))
            node_sequence.append(rnd_untie)

        node_sequence = list(itertools.chain.from_iterable(node_sequence))
        #print(raw_ranking[node_sequence])

    
    return node_sequence
########################################
# functions to attack networks and measure robustness
def erase_node_from_Mnet(Mnet_dict,node_name,**kwargs):


    try:
        by_all=kwargs["by_all"]
    except:
        by_all = False

    verbose=0
    if not (by_all): #only by name
        for node in list(Mnet_dict.keys()):
            if (node[0]==node_name):
                    del Mnet_dict[node]
            else:
                for node2 in list(Mnet_dict[node].keys()):
                    if (node2[0]==node_name):
                        del Mnet_dict[node][node2]
    else:

         for node in list(Mnet_dict.keys()):
            if (node == node_name):
                del Mnet_dict[node]
            else:
                for node2 in list(Mnet_dict[node].keys()):
                    if (node2 == node_name):
                        del Mnet_dict[node][node2]



    if (verbose) :
        print("-- node %s deleted -- " % node_name)

    return Mnet_dict
def measure_secondary_extinctions(Mnet_dict,ext_struct,linking_set):
    verbose=0
    if (verbose) :
        print("entro a medir extinctiones:\n")
        #print(Mnet_dict)
    #first round of secondary extinctions: species directly linked to plants

    for node in list(Mnet_dict.keys()):
        #print("searching for node %s %s %s" % (node[0],node[1],node[2]))
        try:
            if (len(Mnet_dict[node].keys())==0): #if no more neighbours
                if (verbose):
                    print("node %s %s %s extinct" % (node[0],node[1],node[2]))
                interaction=node[2]
                ext_set = node[1]
                #if (verbose):
                    #print(interaction)


                # area_struct["Area_merged"] = \
                #     area_struct["Area_merged"] + \
                #     1./(area_struct["N_passive"]/area_struct["N_active"])

                ext_struct["merged"]= ext_struct["merged"] + 1

                # area_struct["Area_merged_%s" % interaction] = \
                #     area_struct["Area_merged_%s" % interaction] + \
                #     1./(area_struct["N_passive_%s" % interaction]/area_struct["N_active"])
                ext_struct[interaction]= ext_struct[interaction] + 1

                # area_struct["Area_%s" % ext_set] = \
                #     area_struct["Area_%s" % ext_set] + \
                #     1./(area_struct["N_%s" % ext_set]/area_struct["N_active"])
                ext_struct[ext_set] = ext_struct[ext_set]  + 1

                #print("try to erase 2nd nodes")
                if (linking_set != "Plant"):

                    if (interaction == 'herbivory'):
                        node_dual=(node[0],node[1],'parasitism')
                    else :
                        node_dual = (node[0], node[1], 'herbivory')

                    try:
                        Mnet_dict=erase_node_from_Mnet(Mnet_dict,node_dual,by_all=True)

                        if (verbose):
                            print("erased %s %s %s from nterowk" % node_dual[0],node_dual[1],node_dual[2])
                    except:
                        pass
                else:
                    Mnet_dict = erase_node_from_Mnet(Mnet_dict, node, by_all=True)
                #print("resulting in dict:\n")
                #print(Mnet_dict)
                del Mnet_dict[node]
        except:
            pass


    #print("second round!")
    #second round of secondary extinctions: species indirectly linked to plants
    for node in list(Mnet_dict.keys()):
        if (len(Mnet_dict[node].keys())==0): #if no more neighbours
            if (verbose):
                print("2ndary node %s %s %s extinct" % (node[0],node[1],node[2]))
            interaction=node[2]
            if (verbose):
                print(interaction)
            ext_set=node[1]

                # area_struct["Area_merged"] = \
                #     area_struct["Area_merged"] + \
                #     1./(area_struct["N_passive"]/area_struct["N_active"])

            ext_struct["merged"]= ext_struct["merged"] + 1

                # area_struct["Area_merged_%s" % interaction] = \
                #     area_struct["Area_merged_%s" % interaction] + \
                #     1./(area_struct["N_passive_%s" % interaction]/area_struct["N_active"])
            ext_struct[interaction]= ext_struct[interaction] + 1

                # area_struct["Area_%s" % ext_set] = \
                #     area_struct["Area_%s" % ext_set] + \
                #     1./(area_struct["N_%s" % ext_set]/area_struct["N_active"])
            ext_struct[ext_set] = ext_struct[ext_set]  + 1
            del Mnet_dict[node]
    #print("Done obtaining 2nd extinctions -- ")


    if (verbose):
        print(ext_struct)
        print("-------")
        #print(Mnet_dict)
        #print("-------")


    return ext_struct
def calc_ext_area(Mnet,area_struct, ext_struct, node_layers, **kwargs):
    verbose=0
    mode = kwargs["mode"]

    if (mode == "UW"):

        area_struct["Area_merged"] = area_struct["Area_merged"] + \
                                 (ext_struct["merged"]/(area_struct["N_passive"]*area_struct["N_active"]))
    #1./(area_struct["N_passive"]*area_struct["N_active"])

        for my_set in Mnet.slices[1]:
            area_struct["Area_%s" % my_set] = area_struct["Area_%s" % my_set] + \
                                          (ext_struct[my_set]/(area_struct["N_%s" % my_set]*area_struct["N_active"]))
        # 1./(area_struct["N_%s" % ext_set]*area_struct["N_active"])


        for interaction in Mnet.slices[2]:
            area_struct["Area_%s" % interaction] = area_struct["Area_%s" % interaction] + \
                                                      (ext_struct[interaction]/(area_struct["N_passive_%s" %    interaction]*area_struct["N_active"]))
        # 1./(area_struct["N_passive_%s" % interaction]*area_struct["N_active"])

        for interaction in node_layers:
                area_struct["Area_only_%s" % interaction] = area_struct["Area_only_%s" % interaction] + ext_struct[interaction]/(area_struct["N_passive_%s" % interaction]*area_struct["N_active_%s" % interaction])

    if (verbose):
        print(area_struct)
    return area_struct
def calc_ext_area_of_seq(Mnet,nodes_to_erase,active_set,net_name,ranking,nrep,**kwargs):
    verbose=False
    try:
        verbose=kwargs["verbose"]
    except:
        verbose=False
    #print("enter to calc Ext Area")
    #we need to define one ext_are for EACH interaction network/set and another for the whole network. #decision, what happens when we erase a plant that is not connected to species in the set?? I guess you just continue, right? the % of extinct species is for the whole network.
    #Area_struct: N_set: initial number of nodes in each set.
    #             Area_merged: extinction are for the full merged network,
    #             Area_set: extinction area for each set (poock's style)
    #             Area_int: extinction are for each interaction layer separately

    try:
        print_to_file=kwargs["print_to_file"]
    except:
        print_to_file=False

    try:
        return_ext_curve=kwargs["ext_curve"]
    except:
        return_ext_curve=False

    try:
        Area_struct=kwargs["Area_struct"]
        Area_struct=reset_area_struct(Area_struct,Mnet)
    except:
        Area_struct = create_area_struct(Mnet, active_set, nodes_to_erase)

    try:
        Ext_struct=kwargs["Ext_struct"]
        Ext_struct=reset_ext_struct(Ext_struct,Mnet)
    except:
        Ext_struct=create_ext_struct(Mnet)


    try:
        mode=kwargs("mode") # eighet qualitative (all species are equally important) or quantitative (importance of species proportional to abundance)
        if (mode=="cuantitative"):
            active_node_abundances=kwargs["abundances"]

    except:
        mode="qualitative"

    erased_cont=0
    new_Mnet_dict = copy.deepcopy(Mnet._net) #since we destroy the dict we need to copy it each time
    linking_set=get_linking_set(Mnet)
    #if (print_to_file):
    filename="./OUTPUT/Data/Individual_areas/Ext_area_%s_%s_%s.dat" % (net_name,ranking,nrep)

    EA_df=pd.DataFrame(columns=list(Ext_struct))

    for node_name in nodes_to_erase:
        #erase node from network
        node_layer=get_layer_of_node(node_name,Mnet,set_name="Plant") #know the interaction it comes from
        new_Mnet_dict=erase_node_from_Mnet(new_Mnet_dict,node_name)

        if (verbose):
            print("%s node %s (%s) erased, go to measure secondary extinctions:\n" % (erased_cont,node_name,node_layer))

        #measure secondary extinctions
        Ext_struct["Plant"] = erased_cont +1
        for interaction in node_layer:
            Ext_struct["Plant_%s" % interaction] = Ext_struct["Plant_%s" % interaction] + 1
        if (mode == "qualitative"):
            Ext_struct=measure_secondary_extinctions(new_Mnet_dict,Ext_struct,linking_set)
            Area_struct=calc_ext_area(Mnet,Area_struct, Ext_struct, node_layer, mode="UW")

        else:
            Ext_struct=measure_secondary_extinctions(new_Mnet_dict,Ext_struct,linking_set)
            Area_struct = calc_ext_area(Mnet, Area_struct, Ext_struct, mode="W",weigths=active_node_abundances)

        if(print_to_file or return_ext_curve):
            EA_df=EA_df.append(Ext_struct,ignore_index=True)
            #print(Ext_struct)
        erased_cont = erased_cont +1

    EA_df["ext_seq"]=nodes_to_erase
    #print(EA_df)
    #print(EA_df)

    if(print_to_file):
        for i in Area_struct:
            EA_df.loc[:,i]=Area_struct[i]

        EA_df.to_csv(filename)


    if (verbose):
        print("Done ---- ")
        print(Area_struct)
        print(new_Mnet_dict)

    if (return_ext_curve):
        return Area_struct, EA_df
    else:
        return Area_struct