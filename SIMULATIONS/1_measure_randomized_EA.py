from Mnetworks import *
from Metrics import *
from Robustness import *
import sys

network_names=['Sinohara_1_ALL_PH', 'Sinohara_2_ALL_PH', 'Sinohara_3_ALL_PH', 'Sinohara_4_ALL_PH', 'Sinohara_ALL_A_PH', 'Sinohara_ALL_E_PH', 'Sinohara_ALL_I_PH', 'Sinohara_2_E_PH', 'Sinohara_3_E_PH', 'Sinohara_4_I_PH', 'Melian_OO_OO_PH', 'Hackett_1_ALL_PH', 'Hackett_2_ALL_PH', 'Hackett_1_S_PH', 'Hackett_1_GL_PH', 'Pocock_OO_OO_PH',  'Melian_OO_OO_HSD',  'McFayden_ALL_A_HPa', 'McFayden_1_A_HPa', 'McFayden_2_A_HPa', 'McFayden_3_A_HPa', 'McFayden_4_A_HPa', 'McFayden_5_A_HPa', 'McFayden_6_A_HPa', 'McFayden_7_A_HPa', 'McFayden_8_A_HPa', 'McFayden_9_A_HPa', 'McFayden_10_A_HPa', 'McFayden_ALL_B_HPa', 'McFayden_1_B_HPa', 'McFayden_2_B_HPa', 'McFayden_3_B_HPa', 'McFayden_4_B_HPa', 'McFayden_5_B_HPa', 'McFayden_6_B_HPa', 'McFayden_7_B_HPa', 'McFayden_8_B_HPa', 'McFayden_9_B_HPa', 'McFayden_10_B_HPa', 'Hackett_1_ALL_HPa', 'Hackett_1_WL_HPa',   'Melian_OO_OO_PSD','Dattilo_OO_OO_PSD','Dattilo_OO_OO_PA']

Nrep = 1000 + 1  # number of random sequences of specie to delete
Nnet = 100 + 1 #number of randomization of each empirical network
#number of random sequences of specie to delete
#ext_MODE="RND" #RND, OD, ODinv, PR
#null_model="NL" # (NL), (NL2), (K2), (K)
ext_MODE=sys.argv[1]
null_model=sys.argv[2]

for name in network_names:
    Mnet=Read_net_general(name)
    linking_set = get_linking_set(Mnet)
    nodes_to_erase = get_nodes_in_set(Mnet, set_name="Plant")
    #estructures to store extinction area results (1000 for each randomization, and 100 randomization for each empirical network)
    Area_df = Initialize_Area_df(Mnet, Nrep=Nrep, Nnet=Nnet)
    Area_struct = create_area_struct(Mnet, "Plant", nodes_to_erase)
    #Node_impact_df = initialize_node_impact_df(Mnet, "Plant", nodes_to_erase, Area_struct, Nrep=Nrep)
    Ext_struct = create_ext_struct(Mnet)

    for net in range(1, Nnet): #100 randomizations
        print("%s net %s randomization ---------- " % (name, net))
        rnd_Mnet = randomize_pymnet(Mnet, method=null_model)
        # print("randomization OK")
        K_df, new_linking_set = get_full_K_df(name,rnd_Mnet, linking_set=linking_set)

        for rep in range(1, Nrep): #1000 extinction sequences
            species_now = set(K_df.index)
            #in case some nodes are lost in the randomization
            nodes_to_erase_now = list(set(nodes_to_erase).intersection(species_now))
            #get extinction sequence according to extinction scenario
            node_sequence = get_node_sequence(nodes_to_erase_now, K_df.loc[nodes_to_erase_now], ext_MODE)
            #merasure extinction areas of extinction sequence
            ext_area = calc_ext_area_of_seq(rnd_Mnet, node_sequence, linking_set, name, ext_MODE, rep, Area_struct=Area_struct, print_to_file=False)
            #Add EA to dataframe
            Area_df = Add_ext_area_to_Area_df(Area_df, ext_area, rep, net=net)

    # Export extinction area dataframe with the 100X1000 results
    filename = "../OUTPUT/Data/Ext_Area_%s_%s_%s_def.csv" % (name, ext_MODE, null_model)
    Area_df.to_csv(filename)