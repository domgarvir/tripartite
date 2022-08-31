from Mnetworks import *
from Metrics import *
from Robustness import *
import sys
network_names=['Sinohara_1_ALL_PH', 'Sinohara_2_ALL_PH', 'Sinohara_3_ALL_PH', 'Sinohara_4_ALL_PH', 'Sinohara_ALL_A_PH', 'Sinohara_ALL_E_PH', 'Sinohara_ALL_I_PH', 'Sinohara_2_E_PH', 'Sinohara_3_E_PH', 'Sinohara_4_I_PH', 'Melian_OO_OO_PH', 'Hackett_1_ALL_PH', 'Hackett_2_ALL_PH', 'Hackett_1_S_PH', 'Hackett_1_GL_PH', 'Pocock_OO_OO_PH',  'Melian_OO_OO_HSD',  'McFayden_ALL_A_HPa', 'McFayden_1_A_HPa', 'McFayden_2_A_HPa', 'McFayden_3_A_HPa', 'McFayden_4_A_HPa', 'McFayden_5_A_HPa', 'McFayden_6_A_HPa', 'McFayden_7_A_HPa', 'McFayden_8_A_HPa', 'McFayden_9_A_HPa', 'McFayden_10_A_HPa', 'McFayden_ALL_B_HPa', 'McFayden_1_B_HPa', 'McFayden_2_B_HPa', 'McFayden_3_B_HPa', 'McFayden_4_B_HPa', 'McFayden_5_B_HPa', 'McFayden_6_B_HPa', 'McFayden_7_B_HPa', 'McFayden_8_B_HPa', 'McFayden_9_B_HPa', 'McFayden_10_B_HPa', 'Hackett_1_ALL_HPa', 'Hackett_1_WL_HPa',   'Melian_OO_OO_PSD','Dattilo_OO_OO_PSD','Dattilo_OO_OO_PA']

network_names=['Sinohara_2_E_PH']

Nrep=3000 #number of random sequences of specie to delete
 
ext_MODE=sys.argv[1] #RND, DD, ID: Random, decreasing degree, increasing degree

for name in network_names:
    print(name)
    #read network, get linking set and the nodes we are going to erase
    Mnet = Read_net_general(name)
    linking_set=get_linking_set(Mnet)
    nodes_to_erase=get_nodes_in_set(Mnet, set_name="Plant") #we always attack the plants
    K_df, new_linking_set=get_full_K_df(name, Mnet,linking_set=linking_set)
    #create structure to store extinction area etc
    Area_df = Initialize_Area_df(Mnet, Nrep=Nrep)
    Area_struct = create_area_struct(Mnet, "Plant", nodes_to_erase)
    Node_impact_df = initialize_node_impact_df(Mnet, "Plant", nodes_to_erase, Area_struct, Nrep=Nrep)
    Ext_struct = create_ext_struct(Mnet)

    for rep in range (Nrep):
        print("rep %s" % rep, end='\r' )
        # get sequence of extinction depending on the extinction scenario
        node_sequence=get_node_sequence(nodes_to_erase,K_df.loc[nodes_to_erase],ext_MODE)
        #print("calc ext area of seq")
        ext_area, ext_curve =calc_ext_area_of_seq(Mnet, node_sequence, "Plant", name,ext_MODE, rep,Area_struct=Area_struct,ext_curve=True)
        #plot_ext_area_sequence_multiplot(name, ext_MODE, rep)
        #add this particular EA to the dataframe of extinction areas
        Area_df=Add_ext_area_to_Area_df(Area_df,ext_area,rep)
        #and to the node impact df: for obtaining later the importance of plants for the different interaction layers/whole community
        if (ext_MODE == "RND"):
            Node_impact_df=Add_node_seq_to_node_impact_df(Node_impact_df,node_sequence,ext_area,rep)

    # export the extinction area dataframe
    filename = "../OUTPUT/Data/Ext_Area_%s_%s.csv" % (name, ext_MODE)
    Area_df.to_csv(filename)
    # export the node impact dataframe
    if (ext_MODE == "RND"):
        filename2 = "../OUTPUT/Data/Ext_Areas/Node_impact_%s_%s.csv" % (name, ext_MODE)
        Node_impact_df.to_csv(filename2)

