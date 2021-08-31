# -*-coding:Latin-1 -*
import pandas as pd
import numpy as np
import itertools
from pymnet import *
from random import seed, randint, shuffle, sample
from math import ceil
#READING - WRITING MNETWORKS
def Read_net_general(net_name):
    #print(net_name)
    filename="../Networks_DB/%s.csv" % (net_name)
    mnet = MultilayerNetwork(aspects=2)  # group and intereaction type

    with open(filename) as fp:
        first_line = fp.readline()
        for cnt, line in enumerate(fp):
            name1, set1, int_type1, name2 , set2, int_type2, Frequency = line.rstrip().split(",")  # animal, plant, type of interaction
            mnet[name1,set1,int_type1][name2,set2,int_type2]=float(Frequency)
    return mnet

# GET SUBSETS of NODES or LAYERS
#get nodes in the linking set
def get_linking_set(Mnet):

    #this fucntion takes a multinetwork with only two diferent layers of interactions and returns the name of the linking set of nodes:
    if (len(Mnet.slices) >2):
        interactions=list(Mnet.slices[2])
        if (len(interactions)!=2):
            print("ERROR only 2 interactions possible")
            return
        sets_involved = {}
        for interaction in interactions:
            sets_involved[interaction]=[]
        for K1 in Mnet._net.keys():
            sets_involved[K1[2]].append(K1[1])
        for interaction in interactions:
            sets_involved[interaction]=set(sets_involved[interaction])
            #print("interaction:%s - sets=%s" %(interaction,sets_involved[interaction]))
        linking_set=list(sets_involved[interactions[0]].intersection(sets_involved[interactions[1]]))

    #else:
    return linking_set[0]
#get nodes in a specific interaction layer
def get_nodes_in_layer(Mnet,interaction_layer="herbivory",nlinks=False):
    nodes_in_layer=[]
    nLinks=0
    for sp_name, sp_set, sp_layer in Mnet._net.keys():
        if (sp_layer==interaction_layer):
            nodes_in_layer.append(sp_name)
            if (nlinks):
                nLinks += len(Mnet._net[(sp_name,sp_set,sp_layer)])
                #print("%s,%s,%s : %s links" % (sp_name, sp_set, sp_layer,nLinks))
    if (nlinks):
        return nodes_in_layer, nLinks
    else:
        return nodes_in_layer
#get nodes in a specific species set (posibility to also restric to only one interaction layer)
def get_nodes_in_spset(Mnet,sign,species_set="Plant",interaction_layer=False):
    #possible names for species sets are those in the keys of dict_node_simple
    dict_sp_sets=invert_dict(dict_node_simple[sign])
    nodes_in_set=[]
    for sp_name, sp_set, sp_layer in Mnet._net.keys():
        if (interaction_layer):
            if(sp_layer == interaction_layer):
                if sp_set in dict_sp_sets[species_set]:
                    nodes_in_set.append(sp_name)
        else:
            if (sp_set in dict_sp_sets[species_set]):
                nodes_in_set.append(sp_name)

    return nodes_in_set
#get nodes in a specific species set with their full information
def get_full_nodes_in_set(Mnet,set_name="Plant"):
    nodes = []

    for k1 in Mnet._net.keys():
        if (k1[1] == set_name):
            nodes.append(k1)

    return nodes
#get nodes in lsp set without merging smaller sets of species into those defined by the interaction layer
def get_nodes_in_set(Mnet,**kwargs):
#this functionas takes a Mnet in ymnet format and returns a list with all the nodes in a given SET, for all the networks considered in aspect names
    nodes=[]
    set_name=kwargs['set_name']

    #if none aspect is provided we consider all of them
    try :
        aspects=list(Mnet.slices[2]) #this means it has more than one aspect
        if (kwargs['aspect_names']): #only return nodes in of one set in one aspect
            aspects = kwargs['aspect_names']
        else:
            aspects = list(Mnet.slices[2]) #return nodes in one set in different aspects, but we need to keep doubled thow species that are repeated in different aspects, rigth?
        #for each aspectwe need to record all nodes, then erase the repeated and THEN, merge all different aspects (if needed)
        temp_nodes = [[] for i in range(len(aspects))]
        for k1 in Mnet._net.keys():
            for aspect_index in range(len(aspects)):
                if ((k1[1] == set_name) and (k1[2]==aspects[aspect_index])):
                    temp_nodes[aspect_index].append(k1[0])

        #once all nodes sorted, erase repeated from inside each aspec:
        for aspect_index in range(len(aspects)):
            temp_nodes[aspect_index]=list(set(temp_nodes[aspect_index]))

        return list(itertools.chain.from_iterable(temp_nodes))

    except:    #only one aspect!
        for k1 in Mnet._net.keys():
            if (k1[1]==set_name):
                nodes.append(k1[0])

        return list(set(nodes))

    # for aspect in aspects:
    #     for k1 in Mnet._net.keys():
    #         if ((k1[1] == layer_name) and (k1[2]==aspect)):


        # if (len(k1)<3): #not good enough,
        #     if (k1[1]==layer_name):
        #         nodes.append(k1[0])
        # elif (len(k1)==3):
        #     if ((k1[1] == layer_name) and (k1[2]==aspect)):
        #         nodes.append(k1[0])
    return list(set(nodes))

#get the name of the interaction layer one node belongs to
def get_layer_of_node(node_name,Mnet,**kwargs):

    try:
        node_set=kwargs["set_name"]
    except:
        node_set=None
    #print("node is %s, in set %s" % (node_name,node_set))

    if (node_set == "Plant"):
        interactions=[]
        for interaction in Mnet.slices[2]:
            plants_in_layer=get_nodes_in_set(Mnet,set_name=node_set,aspect_names=[interaction])
            #print("plants in %s: %s " % (interaction,plants_in_layer))
            if (node_name in plants_in_layer):
                #print("está aqui!!\n")
                interactions.append(interaction)


    return interactions
#OPERATION IN NETWORKS
# get unilayer network from Mnetwork
def get_1Aspect_network(full_mnt, aspect_name):
    #returns a unimodal (maybe bepartite) network of only one aspect
    sub_net=MultilayerNetwork(aspects=1)

    for k1 in full_mnt._net.keys():
        if (k1[2]==aspect_name):
            for k2 in full_mnt._net[k1].keys():
                sub_net[k1[0],k1[1]][k2[0],k2[1]]=full_mnt[k1[0],k1[1],k1[2]][k2[0],k2[1],k2[2]]

    return sub_net
#add 1 interaction layer to a multinetwork
def add_net_to_Multilayer(Mnet,net_to_add,aspect,**kwargs):

    row_names=False
    try:
        row_layer_name = kwargs['row_layer_name']
        col_layer_name=kwargs['col_layer_name']
        row_names=True
    except:
        pass

    if (row_names):
        for item in net_to_add.keys():
            node_t1=item[0] #lepidoptera
            node_t2=item[1] #plant
            weigth=net_to_add[item]
            #print("link from %s, %s, %s to %s, %s, %s : %s" %(node_t1,row_layer_name,aspect,node_t2,col_layer_name,aspect,weigth))
            Mnet[node_t1,row_layer_name,aspect][node_t2,col_layer_name,aspect]=weigth
    else:
        for nodei in net_to_add._net.keys():
            #nodei_name=nodei[0]
            #nodei_set=nodei[1]
            for nodej in net_to_add._net[nodei].keys():
                #nodej_name=nodej[0]
                #nodej_set=nodej[1]
                weigth=net_to_add[nodei][nodej]
                Mnet[nodei[0],nodei[1],aspect][nodej[0],nodej[1],aspect]=weigth

    return Mnet
#CHANGING NETWORK FORMAT
#merge tripartite network to a bipartite network (linking set species and the rest)
def get_merged_multinetwork(Mnet):

    merged_Mnet = MultilayerNetwork(aspects=1)

    for node in Mnet._net.keys():
        for nnode in Mnet._net[node]:
            try: #if already existent add connections
                merged_Mnet[node[0],node[1]][nnode[0],nnode[1]] = merged_Mnet[node[0],node[1]][nnode[0],nnode[1]] + float(Mnet._net[node][nnode])

            except:

                merged_Mnet[node[0],node[1]][nnode[0],nnode[1]]= float(Mnet._net[node][nnode])

    return merged_Mnet
#merge multinetworks to a new multinetwork with more layers
def merge_multinetworks(Mnets):

    merged_Mnet=MultilayerNetwork(aspects=2)

    for Mnet in Mnets:
        for k1 in Mnet._net.keys():
            for k2 in Mnet._net[k1].keys():
                merged_Mnet[k1[0], k1[1], k1[2]][k2[0], k2[1], k2[2]] = Mnet[k1[0], k1[1], k1[2]][k2[0], k2[1], k2[2]]

    return merged_Mnet
#get multi-incidence matrix form multinetwork
def get_incidence_matrix_from_mnet(net,linking_set="Plant",weigth=False,**kwargs):
    #input: a one layer Mnet
    #for the moment only for one layer, for two is in write_MnetMatrix_new
    #we will have from 2 to 3 indices

    sp_sets=list(net.slices[1])

    net_df=get_df_from_ment(net,weight=weigth) #square version, we need to clean to keep only the incidence part
    #print(net_df.mean().mean())

    if (len(net_df.index.names)==2):
        #print("len 2")
        m1 = net_df.xs(linking_set, level='set', axis=0, drop_level=False)
        sp_sets.remove(linking_set)
        m = m1.T.index.get_level_values('set').isin(sp_sets)
        result = m1.T[m].T

        result = result.reorder_levels(['set', 'name'])
        result = result.reorder_levels(['set', 'name'], axis=1)

    elif ( len(net_df.index.names)==3):
        #print("len 3")
        m1 = net_df.xs(linking_set, level='set', axis=0, drop_level=False)
        m1=m1.groupby(by=['name','set']).sum()
        sp_sets.remove(linking_set)
        m = m1.T.index.get_level_values('set').isin(sp_sets)
        result = m1.T[m].T

        result = result.reorder_levels(['set', 'name'])
        result = result.reorder_levels(['set', 'int','name'], axis=1)

    return result
#get multinet from multi-incidence matrix
def get_net_from_incidence_matrix(net_df,**kwargs):

    try:
        interaction=kwargs["int_layer"]
    except:
        interaction=False

    try:
        multiple_interaction=kwargs["mnet"]
    except:
        multiple_interaction=False

    if (interaction):
        net = MultilayerNetwork(aspects=2)

        row=list(net_df.index)
        col=list(net_df.columns)

        for sp1 in row:
            for sp2 in col:
                net[sp1[1],sp1[0],interaction][sp2[1],sp2[0],interaction]=net_df.loc[sp1,sp2]

    elif (multiple_interaction):
        net = MultilayerNetwork(aspects=2)
        row=list(net_df.index)
        col=list(net_df.columns)

        for sp1 in row:
            for sp2 in col:
                net[sp1[1],sp1[0],sp2[1]][sp2[2],sp2[0],sp2[1]]=net_df.loc[sp1,sp2]

    return net
#change mnetwork from pymnet format to multiple nx unilayer networks
def from_pymnet_2_nx(Mnet,**kwargs):
#this functions takes a multilayer network in the pymnet format and returns a dictionary with as many keys as networks in the original mnet, each contanining the "aspect", and in each value the unimodal network in nx format.
    aspects=Mnet.slices[2] #this are all the unimodal networks contained in the Mnet
    nxMnet={}
    for aspect in aspects:#for each network/aspect
        bipartite_sets = []  # these are the different sets in the network. In an unimodal we will only find 2 of them
        #print("aspect:%s" %aspect)
        temp_net=nx.Graph()
        for k1 in Mnet._net.keys():#for each node
            if (k1[2] == aspect): #if in the chosen network/aspect
                for k2 in Mnet._net[k1].keys():#for each neighbour
                    temp_net.add_edge(k1[0],k2[0],weight=Mnet[k1[0],k1[1],k1[2]][k2[0],k2[1],k2[2]])
                temp_net.nodes[k1[0]]['type']=k1[1] #add name of set the node is in
                #print("node:%s type:%s from BS:%s" %(k1[0],k1[1],bipartite_sets))
                if (k1[1] not in bipartite_sets):
                    bipartite_sets.append(k1[1])
                if (k1[1] == bipartite_sets[0]):
                    temp_net.nodes[k1[0]]['bipartite'] =0
                else:
                    temp_net.nodes[k1[0]]['bipartite'] =1
                nxMnet[aspect]=temp_net
        #print("aspect: %s  bipertite sets:%s" %(aspect,bipartite_sets))
    return nxMnet
#and reverse change
def from_nx_2_pymnet(nxMnet,**kwargs):
#this functions takes a dictinary contaning a multilayer network, each network in one value, and the key containing the "aspect" name (type of network), and returns the multilayer network in the pymnet format.
    aspects=list(nxMnet.keys())
    mnet = MultilayerNetwork(aspects=2)  # multipartite, contains information about interaction and partition set.

    for aspect in aspects:
        edges=list(nxMnet[aspect].edges())
        for edge in edges:
            mnet[edge[0],nxMnet[aspect].nodes[edge[0]]['type'],aspect][edge[1],nxMnet[aspect].nodes[edge[1]]['type'],aspect]=nxMnet[aspect].edges[edge]['weight']

    return mnet
#from pymnet to matrix
def get_df_from_ment(Mnet,**kwargs):

    weigth=False
    try:
        weigth=kwargs["weight"]
    except:
        pass

    multipartite=True
    try:
        sets=Mnet.slices[1]
        interactions=Mnet.slices[2]
        indices=['name','set','int','k']
    except:
        print("no multipartite")
        sets = Mnet.slices[1]
        indices = ['name', 'set', 'k']
        multipartite=False

    nodes = list(Mnet._net.keys())
    nodes_with_k=[]
    for node in nodes:
        if (multipartite):
            nodes_with_k.append((node[0],node[1],node[2],len(Mnet._net[node].keys())))
            mindex = pd.MultiIndex.from_tuples(nodes_with_k, names=('name', 'set', 'int', 'k'))
        else:
            nodes_with_k.append((node[0], node[1], len(Mnet._net[node].keys())))
            mindex = pd.MultiIndex.from_tuples(nodes_with_k, names=('name', 'set', 'k'))


    # mindex=mindex.sortlevel(level = ['int','set'])
    M_df=pd.DataFrame(0, index=mindex, columns=mindex)
    #M_df.to_csv("caca.csv")
    for nodei in Mnet._net.keys():
        #print(nodei[0], nodei[1], nodei[2])
        for nodej in Mnet._net[nodei].keys():
            if (multipartite):
                keyi = (nodei[0], nodei[1], nodei[2], len(Mnet._net[nodei].keys()))
                keyj = (nodej[0], nodej[1], nodej[2], len(Mnet._net[nodej].keys()))
            else :
                keyi = (nodei[0], nodei[1], len(Mnet._net[nodei].keys()))
                keyj = (nodej[0], nodej[1], len(Mnet._net[nodej].keys()))

            if (weigth):
                #M_df.loc[nodei,nodej]=Mnet._net[nodei][nodej]
                #M_df.loc[nodej,nodei]=Mnet._net[nodej][nodei]
                M_df.loc[keyi, keyj] = float(Mnet._net[nodei][nodej])*1.
                M_df.loc[keyj, keyi] = float(Mnet._net[nodej][nodei])*1.
            else:
                M_df.loc[keyi, keyj] = 1
                M_df.loc[keyj, keyi] = 1

    #if(weigth):


    if multipartite:
        M_df.sort_index(axis=0,level=['set','int','k'], ascending = [False, False, False], inplace=True)
        M_df.sort_index(axis=1, level=['set', 'int','k'], ascending = [False, False, False], inplace=True)

    else :
        M_df.sort_index(axis=0, level=['set', 'k'], ascending=[False, False], inplace=True)
        M_df.sort_index(axis=1, level=['set', 'k'], ascending=[False, False], inplace=True)

    M_df = M_df.droplevel('k', axis=1).droplevel('k', axis=0)

    #print(M_df.head(15))

    return M_df
#NETWORK RANDOMIZATIONS
############################################
def nullmodel_constantNL(nxMnet, **kwargs):
    # this null model only maintains the same number of links in each network, irrespective of the node degree. If nothing is said, all aspects are randomized in te same fashion
    aspects = list(nxMnet.keys())#interactions
    try:
        aspects = kwargs['aspects_to_rewire']  # we can select the aspects to rewire
    except:
        pass
    rnd_mnet = {} #returns a nxMnet format too

    for aspect in aspects:
        #print("aspect %s:" % aspect)
        rnd_net=rewire_constantNL(nxMnet[aspect])
        rnd_mnet[aspect]=rnd_net

    return rnd_mnet
def rewire_constantNL(net):
    rnd_net = nx.Graph()
    #first copy all original nodes so they remain there:
    rnd_net.add_nodes_from(net.nodes(data=True))
    #in this rewrirign we will create as many links in the new networks as where in the original, irrespective of the degree:
    change=0
    nodes_a=[]
    nodes_b=[]
    #get different sets of nodes, since network may be not completely connected we sort like this, based in the label
    for node in net.nodes():
        if (net.nodes[node]['bipartite']):
            nodes_b.append(node)
        else:
            nodes_a.append(node)

    edges=list(rnd_net.edges())

    #copy lists because we are going to force that each node has at least 1 connection
    unconnected_nodes_a = nodes_a.copy()
    unconnected_nodes_b = nodes_b.copy()

    #print("comenzamos con %s edges" % len(list(net.edges())))
    while (change < len(list(net.edges()))):


        #print("change %s of %s" %(change,len(list(net.edges()))))
        from_unconnected_a=False
        if (len(unconnected_nodes_a)>0): #while there are unconnected nodes
            node_a=unconnected_nodes_a[randint(0,len(unconnected_nodes_a)-1)]
            from_unconnected_a = True
        else:
            node_a = nodes_a[randint(0, len(nodes_a) - 1)]

        from_unconnected_b=False
        if (len(unconnected_nodes_b)>0):
            node_b = unconnected_nodes_b[randint(0, len(unconnected_nodes_b) - 1)]
            from_unconnected_b=True
        else:
            node_b=nodes_b[randint(0,len(nodes_b)-1)]

        while ( ((node_a,node_b) in edges) or ((node_b,node_a) in edges) ):
            #in principle it is not possible to have repeated edges if taken from unconnected nodes
            node_b=nodes_b[randint(0,len(nodes_b)-1)]
            node_a = nodes_a[randint(0,len(nodes_a) - 1)]

        if (from_unconnected_a):
            unconnected_nodes_a.remove(node_a)
        if (from_unconnected_b):
            unconnected_nodes_b.remove(node_b)

        rnd_net.add_edge(node_a,node_b,weight=1)#unweighted
        edges = list(rnd_net.edges())
        #print(edges)
        change = change +1

    return rnd_net
def nullmodel_constantK(nxMnet,**kwargs):
    # this null model maintains the node degree in each network, but changes the neigbour it is linked to. If nothis is said all aspects are randomized in the same fashion
    aspects = list(nxMnet.keys())
    try:
        aspects = kwargs['aspects_to_rewire']  # we can select the aspects to rewire
    except:
        pass

    rnd_mnet = {}  # returns a nxMnet format too

    for aspect in aspects:
        #print("rewire %s network:" % (aspect))
        rnd_net = rewire_constantK(nxMnet[aspect])
        rnd_mnet[aspect] = rnd_net

    return rnd_mnet
def rewire_constantK(net): #unweighted
    rnd_net = net.copy()
    edges=list(rnd_net.edges())
    #print("comenzamos con %s edges:"%len(edges))
    #print(edges)
    #we make as many changes as links:
    change=0
    while (change<len(list(net.edges()))-1): #len(list(net.edges()))-1
        #print("change %s of %s" %(change,len(edges)))
        #print(edges)
        edge1=edges[randint(0,len(edges)-1)]
        edge2=edges[randint(0,len(edges)-1)]
        #print("edge1=%s edge2=%s" % (edge1, edge2))
        while ((edge2 == edge1) and (len(edges)>3)): #we make sure we dont select the same edge twice
            #print("change edge2!")
            edge2=edges[randint(0,len(edges)-1)]
        #print("final edge1=%s edge2=%s" % (edge1, edge2))
        n1a,n1b=edge1
        n2a,n2b=edge2
        #since we are working with bipartite we need to ensuer that links are between different sets!
        if (rnd_net.nodes[n1a]['type'] != rnd_net.nodes[n2a]['type']):
            edge2=n2b,n2a
            #print("change order!")
        #print("n1a=%s n1b=%s n2a=%s n2b=%s" %(n1a,n1b,edge2[0],edge2[1]))
        #we look if the new edges already exists:
        #print("looking for edge %s-%s and %s-%s" %(n1a,edge2[1],edge2[0],n1b))
        #print("%s & %s AND %s & %s =%s " % (((n1a,edge2[1]) not in edges), ((edge2[1],n1a) not in edges),((edge2[0],n1b) not in edges),((n1b,edge2[0]) not in edges), (((n1a,edge2[1]) not in edges) and ((edge2[1],n1a) not in edges)) and (((edge2[0],n1b) not in edges) and ((n1b,edge2[0]) not in edges))   ))
        if (  (((n1a,edge2[1]) not in edges) and ((edge2[1],n1a) not in edges)) and (((edge2[0],n1b) not in edges) and ((n1b,edge2[0]) not in edges))   ) :
            #print(list(rnd_net.edges()))
            rnd_net.remove_edge(n1a,n1b)
            rnd_net.remove_edge(n2a,n2b)
            #print("removing edges %s and %s" %(edge1,edge2))
            #print(list(rnd_net.edges()))
            rnd_net.add_edge(n1a,edge2[1],weight=1)
            rnd_net.add_edge(edge2[0],n1b,weight=1)
            #print("added edges %s-%s and %s-%s" % (n1a,edge2[1],edge2[0],n1b))
            edges = list(rnd_net.edges()) #actualize the edge list
            #print("changed!")
            #print(edges)
            change = change +1
        #else:
            #print("repeated, do nothing")
        #change = change + 1

    #print("terminamos con %s edges:" % len(list(rnd_net.edges())))
    #print(rnd_net.edges())
    return rnd_net
def find_presences(input_matrix):
    num_rows, num_cols = input_matrix.shape
    hp = []
    iters = num_rows if num_cols >= num_rows else num_cols
    input_matrix_b = input_matrix if num_cols >= num_rows else np.transpose(input_matrix)
    for r in range(iters):
        hp.append(list(np.where(input_matrix_b[r] == 1)[0]))

    return hp
def curve_ball(input_matrix, r_hp, num_iterations=-1):
    num_rows, num_cols = input_matrix.shape
    l = range(len(r_hp))
    num_iters = 5*min(num_rows, num_cols) if num_iterations == -1 else num_iterations
    for rep in range(num_iters):
        AB = sample(l, 2)
        a = AB[0]
        b = AB[1]
        ab = set(r_hp[a])&set(r_hp[b]) # common elements
        l_ab=len(ab)
        l_a=len(r_hp[a])
        l_b=len(r_hp[b])
        if l_ab not in [l_a,l_b]:
            tot=list(set(r_hp[a]+r_hp[b])-ab)
            ab=list(ab)
            shuffle(tot)
            L=l_a-l_ab
            r_hp[a] = ab+tot[:L]
            r_hp[b] = ab+tot[L:]
    out_mat = np.zeros(input_matrix.shape, dtype='int8') if num_cols >= num_rows else  np.zeros(input_matrix.T.shape, dtype='int8')
    for r in range(min(num_rows, num_cols)):
        out_mat[r, r_hp[r]] = 1
    result = out_mat if num_cols >= num_rows else out_mat.T
    return result
def nullmodel_constantKn(mnet,**kwargs):
    #print("new type of random model")
    linking_set = get_linking_set(mnet)
    interactions = list(mnet.slices[2])
    rnd_Mnet=MultilayerNetwork(aspects=2)
    for interaction in interactions:
        #print(interaction)
        net = get_1Aspect_network(mnet, interaction)
        net_df=get_incidence_matrix_from_mnet(net, linking_set=linking_set)
        r_hp=find_presences(net_df.values)
        RM=curve_ball(net_df.values,r_hp) #randomized by curveball
        rnd_net_df=pd.DataFrame(RM,index=net_df.index, columns=net_df.columns)
        rnd_net=get_net_from_incidence_matrix(rnd_net_df,int_layer=interaction)
        rnd_Mnet = merge_multinetworks([rnd_net,rnd_Mnet])
    return rnd_Mnet

# def nullmodel_constantKnonLS_constantN(mnet):
#     #this functions randomizes the two layers separately, keeping constant the ki of the non LS nodes, and preventing LSnodes that are not in one interaction layer entering it.
#     linking_set=get_linking_set(mnet)
#     rnd_mnet= MultilayerNetwork(aspects=2)
#     interactions = list(mnet.slices[2])
#     #randomixe independently each interaction layer
#     for interaction in interactions:
#         sub_mnet=get_1Aspect_network(mnet,interaction)
#         sub_df=get_incidence_matrix_from_mnet(sub_mnet,linking_set=linking_set)
#         #sub_df_T=sub_df.T
#         #randomize the rows (keep constant ki  of nonLS nodes)
#         #np.random.shuffle(sub_df.values)
#         #we need to shuffle each column independently
#         is_disconnected=True
#         while is_disconnected:
#             rnd_df=sub_df.copy(deep=True)
#             for col in range(rnd_df.shape[1]):
#                 np.random.shuffle(rnd_df.T.values[col])
#
#             if (not 0 in rnd_df.sum(axis=1).values):
#                 is_disconnected=False
#
#         rnd_sub_mnet=get_net_from_incidence_matrix(rnd_df,int_layer=interaction)
#         #now merge both to return the full randomized network
#         rnd_mnet = add_net_to_Multilayer(rnd_mnet, rnd_sub_mnet, interaction)
#
#     return rnd_mnet
def nullmodel_constantKnonLS_constantN(nxMnet, linking_set,**kwargs):
    aspects = list(nxMnet.keys())  # interactions
    try:
        aspects = kwargs['aspects_to_rewire']  # we can select the aspects to rewire
    except:
        pass
    rnd_mnet = {}  # returns a nxMnet format too

    for aspect in aspects:
        #print("aspect %s:" % aspect)
        rnd_net = rewire_constantKnonLS_constantNL(nxMnet[aspect],linking_set)
        rnd_mnet[aspect] = rnd_net

    return rnd_mnet
def rewire_constantKnonLS_constantNL(net,linking_set):
    #print("enter to rewire, linking set: %s" % linking_set)
    rnd_net = nx.Graph()
    rnd_net.add_nodes_from(net.nodes(data=True))

    nodes_a = []
    unconnected_nodes_a=[]
    nodes_b = []
    unconnected_nodes_b=[]
    dictionary_of_nonLS_degree={}

    # get different sets of nodes, since network may be not completely connected we sort like this, based in the label
    for node in net.nodes():
        if (net.nodes[node]['bipartite']):
            if (net.nodes[node]['type'] == linking_set): #nodes in linking set only included once
                which_linking_set = "nodes_b"
                which_no_linking_set="nodes_a"
                nodes_b.append(node)
            else:
                dictionary_of_nonLS_degree[node] = net.degree(node)
                #print("node %s has %s links so we add it %s times" % (node,net.degree(node),net.degree(node)-1))
                for i in range(net.degree(node)-1): #if not in linking set add as many times as links
                    #print("append!")
                    nodes_b.append(node)

            unconnected_nodes_b.append(node)
        else:
            if(net.nodes[node]['type'] == linking_set):
                which_linking_set = "nodes_a"
                which_no_linking_set = "nodes_b"
                nodes_a.append(node)
            else:
                dictionary_of_nonLS_degree[node] = net.degree(node)
                #print("node %s has %s links so we add it %s times" % (node, net.degree(node), net.degree(node) - 1))
                for i in range(net.degree(node)-1): #if not in linking set add as many times as links
                    #print("append!")
                    nodes_a.append(node)
            unconnected_nodes_a.append(node)

    unconnected_nodes={"nodes_a":unconnected_nodes_a,"nodes_b":unconnected_nodes_b}
    nodes={"nodes_a":nodes_a,"nodes_b":nodes_b}

   #we start by the nodes in the linking set: choose at least one random partner for each one:
    Nlinks=len(net.edges())
    available_links=Nlinks
    #print("we need to assing %s links" % Nlinks)
    #print("we have %s and %s links in nonLSspecies" % (len(unconnected_nodes[which_no_linking_set]),len(nodes[which_no_linking_set])))

    #for the frist round of assignations we dont need to check anything
    from_unconnected_nodes = True
    for node_ls in unconnected_nodes[which_linking_set]:
        nb_of_unconneted=len(unconnected_nodes[which_no_linking_set])

        if (nb_of_unconneted>0):
            node_nls=unconnected_nodes[which_no_linking_set][randint(0,nb_of_unconneted -1)]
        else:
            nb_of_possible_partners=len(nodes[which_no_linking_set])
            node_nls=nodes[which_no_linking_set][randint(0,nb_of_possible_partners -1)]
            from_unconnected_nodes = False

        rnd_net.add_edge(node_ls, node_nls, weight=1)  # unweighted
        if (from_unconnected_nodes):
            unconnected_nodes[which_no_linking_set].remove(node_nls)
            #unconnected_nodes_a.remove(node_a)
        else:
            nodes[which_no_linking_set].remove(node_nls)

        available_links -= 1
    #print("after first round of LS assignation still %s links remaining" % available_links)
    #so now we have all linking set species connected at least one, now we need to connect the rest of non linking_set species with the remaining_links
    #print("still remains %s unconnected NonLS nodes" % len(unconnected_nodes[which_no_linking_set]))
    for node_nls in unconnected_nodes[which_no_linking_set]:
        nb_of_possible_partners = len(nodes[which_linking_set])
        node_ls=nodes[which_linking_set][randint(0,nb_of_possible_partners -1)]

        rnd_net.add_edge(node_ls, node_nls, weight=1)  # unweighted
        available_links -= 1
    #print("after second round of nonLS assignation still %s links remaining" % available_links)
    #now we have everyone connected, so now only remains to distribute  the available links according to the degree of the non-ls species:
    #remaining_species=len(nodes[which_no_linking_set]) #this should be the same as the nb of available links:
    remaining_dict = {i:nodes[which_no_linking_set].count(i) for i in nodes[which_no_linking_set]}
    #print(remaining_dict)
    for node_nls, links in remaining_dict.items():
        available_species=nodes[which_linking_set].copy()
        neighbors=nx.neighbors(rnd_net, node_nls)
        for nn in neighbors:
            available_species.remove(nn)
        #now select the rest of partners freely:
        for i in range(links):
            node_ls=available_species[randint(0,len(available_species)-1)]
            rnd_net.add_edge(node_ls, node_nls, weight=1)  # unweighted
            available_species.remove(node_ls)

    return rnd_net

def nullmodel_constantKnonLS_constantPkLS(mnet):
    linking_set=get_linking_set(mnet)
    rnd_mnet = MultilayerNetwork(aspects=2)
    interactions = list(mnet.slices[2])
    for interaction in interactions:
        sub_mnet=get_1Aspect_network(mnet,interaction)
        sub_df=get_incidence_matrix_from_mnet(sub_mnet,linking_set=linking_set)
        #now we need to shuffle the values of the index_level "name"
        sp_names=list(sub_df.index.get_level_values("name"))
        rnd_sp_names=np.random.permutation(sp_names)
        rename_dict = dict(zip(sp_names, rnd_sp_names))
        #randomly change names of LS nodes
        sub_df.rename(index=rename_dict,level="name",inplace=True)
        rnd_sub_mnet = get_net_from_incidence_matrix(sub_df, int_layer=interaction)
        rnd_mnet = add_net_to_Multilayer(rnd_mnet, rnd_sub_mnet, interaction)

    return rnd_mnet

def randomize_pymnet(Mnet,**kwargs):

    method=kwargs["method"] #possibilities: L (constant number of links) , K (constant connectivity)

    nxMnet = from_pymnet_2_nx(Mnet)
    methods_with_nx=["NL","constant_L","constant_K","NL2"]
    #print("get nxMnet")
    if (method=="NL"):
        rnd_nxMnet = nullmodel_constantNL(nxMnet)
    elif (method=="constant_K"):
        rnd_nxMnet = nullmodel_constantK(nxMnet)
    elif (method=="K"):
        rnd_Mnet = nullmodel_constantKn(Mnet)
    elif (method=="NL2"):
        #rnd_Mnet = nullmodel_constantKnonLS_constantN(Mnet)
        linking_set = get_linking_set(Mnet)
        rnd_nxMnet = nullmodel_constantKnonLS_constantN(nxMnet, linking_set)
    elif (method) == ("K2"):
        rnd_Mnet = nullmodel_constantKnonLS_constantPkLS(Mnet)
    else:
        print("ERROR in null model")

    if (method in methods_with_nx):
        rnd_Mnet = from_nx_2_pymnet(rnd_nxMnet)
        #print("go back to rnd Mnet")

    return rnd_Mnet

#Helpers
#from interaction type to sign of network
dict_net_type={"H-P":"MA", "H-SD":"MA", "H-Pa":"AA", "P-SD":"MM","P-A":"MM"}
#from net_name to interaction type
dict_name_net={
                #MA
               "Sinohara_1_ALL_PH":"H-P","Sinohara_2_ALL_PH":"H-P", "Sinohara_3_ALL_PH":"H-P","Sinohara_4_ALL_PH":"H-P",
               "Sinohara_ALL_A_PH":"H-P", "Sinohara_ALL_E_PH":"H-P", "Sinohara_ALL_I_PH":"H-P",
               "Sinohara_1_A_PH":"H-P","Sinohara_2_A_PH":"H-P", "Sinohara_3_A_PH":"H-P","Sinohara_4_A_PH":"H-P",
               "Sinohara_1_E_PH":"H-P","Sinohara_2_E_PH":"H-P", "Sinohara_3_E_PH":"H-P","Sinohara_4_E_PH":"H-P",
               "Sinohara_1_I_PH":"H-P","Sinohara_2_I_PH":"H-P", "Sinohara_3_I_PH":"H-P","Sinohara_4_I_PH":"H-P",
               "Melian_OO_OO_PH":"H-P",
               "Hackett_1_ALL_PH":"H-P", "Hackett_2_ALL_PH":"H-P",
               "Hackett_1_SM_PH":"H-P_old","Hackett_1_S_PH":"H-P","Hackett_1_SD_PH":"H-P_old","Hackett_1_WL_PH":"H-P_old",
               "Hackett_1_GL_PH":"H-P","Hackett_1_HL_PH":"H-P_old",
               "Hackett_2_SD_PH":"H-P_old",  "Hackett_2_SM_PH":"H-P_old", "Hackett_2_SC_PH":"H-P_old",
               "Melian_OO_OO_HSD":"H-SD",
               "Pocock_OO_OO_PH":"H-P",
                #AA
               "McFayden_1_A_HPa":"H-Pa", "McFayden_2_A_HPa":"H-Pa", "McFayden_3_A_HPa":"H-Pa", "McFayden_4_A_HPa":"H-Pa", "McFayden_5_A_HPa":"H-Pa", "McFayden_6_A_HPa":"H-Pa", "McFayden_7_A_HPa":"H-Pa", "McFayden_8_A_HPa":"H-Pa", "McFayden_9_A_HPa":"H-Pa", "McFayden_10_A_HPa":"H-Pa", "McFayden_ALL_B_HPa":"H-Pa", "McFayden_1_B_HPa":"H-Pa", "McFayden_2_B_HPa":"H-Pa", "McFayden_3_B_HPa":"H-Pa", "McFayden_4_B_HPa":"H-Pa", "McFayden_5_B_HPa":"H-Pa", "McFayden_6_B_HPa":"H-Pa", "McFayden_7_B_HPa":"H-Pa", "McFayden_8_B_HPa":"H-Pa", "McFayden_9_B_HPa":"H-Pa", "McFayden_10_B_HPa":"H-Pa",
               "Hackett_1_ALL_HPa":"H-Pa", "Hackett_2_ALL_HPa":"H-Pa_old",
"McFayden_1_A_PH":"H-Pa", "McFayden_2_A_PH":"H-Pa", "McFayden_3_A_PH":"H-Pa", "McFayden_4_A_PH":"H-Pa", "McFayden_5_A_PH":"H-Pa", "McFayden_6_A_PH":"H-Pa", "McFayden_7_A_PH":"H-Pa", "McFayden_8_A_PH":"H-Pa", "McFayden_9_A_PH":"H-Pa", "McFayden_10_A_PH":"H-Pa", "McFayden_ALL_B_PH":"H-Pa", "McFayden_1_B_PH":"H-Pa", "McFayden_2_B_PH":"H-Pa", "McFayden_3_B_PH":"H-Pa", "McFayden_4_B_PH":"H-Pa", "McFayden_5_B_PH":"H-Pa", "McFayden_6_B_PH":"H-Pa", "McFayden_7_B_PH":"H-Pa", "McFayden_8_B_PH":"H-Pa", "McFayden_9_B_PH":"H-Pa", "McFayden_10_B_PH":"H-Pa",
               "Hackett_1_ALL_PH":"H-P", "Hackett_2_ALL_PH":"H-P",
    "Hackett_1_SM_HPa":"H-Pa", #"Hackett_1_SD_HPa":"H-Pa", "Hackett_2_SD_HPa":"H-Pa",
    "Hackett_1_S_HPa":"H-Pa_old",
    "Hackett_1_GL_HPa":"H-Pa_old","Hackett_1_HL_HPa":"H-Pa_old", "Hackett_1_WL_HPa":"H-Pa",
               "Hackett_2_SC_HPa":"H-Pa_old", "Hackett_1_RB_HPa":"H-Pa_old",
               "Hackett_1_ALL_SHPa":"SH-Pa", "Hackett_2_ALL_SHPa":"SH-Pa",
    "Hackett_1_SM_HPa":"H-Pa", "Hackett_1_WL_SHPa":"SH-Pa",
               "McFayden_ALL_A_HPa":"H-Pa", "McFayden_ALL_A_PH":"H-Pa",
               "Hackett_1_S_SHPa":"SH-Pa","Hackett_1_SD_SHPa":"SH-Pa",  "Hackett_2_SD_SHPa":"SH-Pa",
               "Hackett_1_GL_SHPa":"SH-Pa","Hackett_1_HL_SHPa":"SH-Pa",  "Hackett_2_SC_SHPa":"SH-Pa",
    "Hackett_1_RB_SHPa":"SH-Pa",
               "Hackett_1_ALL_LHSH":"LH-SH","Hackett_2_ALL_LHSH":"LH-SH", "Hackett_1_S_LHSH":"LH-SH","Hackett_1_SD_LHSH":"LH-SH", "Hackett_2_SD_LHSH":"LH-SH","Hackett_1_GL_LHSH":"LH-SH","Hackett_1_HL_LHSH":"LH-SH", "Hackett_2_SC_LHSH":"LH-SH",
                #MM
               "Melian_OO_OO_PSD":"P-SD",
               "Dattilo_OO_OO_PSD":"P-SD", "Dattilo_OO_OO_SDA":"SD-A", "Dattilo_OO_OO_PA":"P-A"
                }
#list of `positive' interactions
positive_interactions=["dispersion", "pollination","nectarivory","Mutualistic","ant","J"]
#list of negative interactions
negative_interactions=["seed_herbivory", "leaf_herbivory","herbivory","frugivory","Trophic","parasitism","H"]
#dictionaries of original species sets names to uniform species set names
dict_node_simple_set_AA={
    'leaf_miner':"Host",
    'leaf_miner_parasitoid':'Parasitoid',
    "Plant":"Plant",
    "Host": "Host",
    "Parasitoid":"Parasitoid"}
dict_node_simple_set_MA={
    'Insect_h':"Herbivore",
    'Insect_p': "Pollinator",
    'flower_visitor':"Pollinator",
    'leaf_miner':"Herbivore",
    'seed_feeder':"Herbivore",
    "flower.visitor":"Pollinator",
    "butterfly":"Pollinator",
    "Disperser":"Seed disperser",
    "Floral visitor":"Pollinator",
    "seed-feeding.rodent":"Herbivore",
    "seed-feeding.insect":"Herbivore",
    "seed-feeding.bird":"Herbivore",
    "aphid":"Herbivore",
    "Plant":"Plant",
    "Herbivore": "Herbivore",
    "Pollinator": "Pollinator"}
dict_node_simple_set_MM={
    'Disperser':"Seed disperser",
    'Seed dispersal':"Seed disperser",
    'flower_visitor':"Pollinator",
    'leaf_miner':"Herbivore",
    'seed_feeder':"Herbivore",
    "flower.visitor":"Pollinator",
    "butterfly":"Pollinator",
    "Floral visitor":"Pollinator",
    "seed-feeding.rodent":"Herbivore",
    "seed-feeding.insect":"Herbivore",
    "seed-feeding.bird":"Herbivore",
    "aphid":"Herbivore",
    "Plant":"Plant",
    "Pollinator":"Pollinator",
    "Herbivore": "Herbivore",
    "Ant":"Ant"}
dict_node_simple={"MA":dict_node_simple_set_MA,"AA":dict_node_simple_set_AA,"MM":dict_node_simple_set_MM}
#dictionaries to obtain the sp name set from the interaction layer name
dict_sp_set_name_from_interaction_M={"pollination":"Pollinator","herbivory":"Herbivore",'dispersion':"Seed disperser","ant":"Ant"}
dict_sp_set_name_from_interaction_A={'parasitism':"Parasitoid",'herbivory':"Plant",'dispersion':"Seed disperser"}
dict_sp_set_name_from_interaction={"MA": dict_sp_set_name_from_interaction_M,"MM": dict_sp_set_name_from_interaction_M, "AA":dict_sp_set_name_from_interaction_A}

cite_dict={"Melian":"\\cite{Melian2009}", "Sinohara": "\cite{Shinohara2019}", "Pocock": "\cite{Pocock2012}","Hackett": "\cite{Hackett2019}", "McFayden":"\cite{Macfadyen2009}", "Dattilo": "\cite{Dattilo2016}"}

#invert dictionary function
def invert_dict(my_dict):
    inv_map = {v:[k for k in my_dict if my_dict[k] == v] for v in my_dict.values()}
    return inv_map
#Invert nested dictionary
def invert_dict_wl(my_dict):
    new_dict = {}
    for k, v in my_dict.items():
        if (isinstance(v,str)):
            new_dict[v] = k
        else:
            for element in v:
                new_dict[element]=k

    return new_dict

