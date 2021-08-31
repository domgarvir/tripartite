from Mnetworks import *
from Metrics import *

import matplotlib.pyplot as plt
import string
import seaborn as sns
from statannot import add_stat_annotation # to annotate boxplots
from scipy import stats

def plot_boxplot_fromdf_wseaborn_to_ax(emp_df,metric,ax=None,**kwargs):
    #color dict for boxplot fill as a function of the sign of networks or their network type
    ax = ax or plt.gca()

    try:
        box_by=kwargs["box_by"]
    except:
        box_by="SIGN" ##SIGN

    if (box_by == "INT"):
        my_column="int"
        columns=["H-P","H-SD","H-Pa","P-SD","P-A"]
        box_pairs =[("H-Pa","H-P"),("P-SD","H-P"),("P-SD","H-Pa")]

    elif (box_by == "SIGN"):
        my_column = "sign"
        columns=["AA","MA","MM"]
        box_pairs=[("AA", "MA") , ("AA", "MM"),("MA", "MM")]

    else:
        my_column=box_by
        columns=emp_df[box_by].unique()
        print(columns)
        if (box_by=="int"):
            box_pairs =[('herbivory','pollination'),('herbivory','dispersion'),('herbivory','parasitism'),('pollination','dispersion'), ('pollination','parasitism') ,('dispersion', 'parasitism')]
        elif (box_by == "set"):
            box_pairs =[('Pollinator', 'Plant'), ('Pollinator', 'Seed disperser'), ('Pollinator', 'Herbivore'), ('Pollinator', 'Host'), ('Pollinator', 'Parasitoid'),   ('Plant', 'Seed disperser'), ('Plant', 'Herbivore'), ('Plant', 'Host'), ('Plant', 'Parasitoid'),   ('Seed disperser', 'Herbivore'), ('Seed disperser', 'Host'), ('Seed disperser', 'Parasitoid'),  ('Herbivore', 'Host'), ('Herbivore', 'Parasitoid'), ('Host', 'Parasitoid')]

        else:
            box_pairs=[("AA", "MA") , ("AA", "MM"),("MA", "MM")]


    try:
        ext_mode=kwargs["ext_mode"]
    except:
        ext_mode=None

    color=True
    try:
        types_color = [color_dict[x] for x in columns]
    except:
        color=False

    annot=False
    try:
        annot=kwargs["annot"]
    except:
        pass

    points=False
    try:
        points=kwargs["points"]
    except:
        pass

    significance_lines=False
    try:
        significance_lines=kwargs["significance_lines"]
    except:
        pass

    null_model = False
    try:
        null_model=kwargs["null_model"]
    except:
        pass

    index = None
    try:
        index=kwargs["index"]
    except:
        pass

    title = None
    try:
        title=kwargs["title"]
    except:
        pass

    ylabel = None
    try:
        ylabel=kwargs["ylabel"]
    except:
        pass

    y_lim = None
    try:
        y_lim=kwargs["y_lim"]
    except:
        pass

    fontsize = 21
    try:
        fontsize = kwargs["fontsize"]
    except:
        pass
    plt.rcParams.update({'font.size': fontsize})

    sns.set(font_scale=1.2)
    bplot = sns.boxplot(data=emp_df, x=my_column, y=metric, width=0.5,order=columns, ax=ax,sym='')
    # change color
    if (color):
        for i in range(0, len(columns)):
            mybox = bplot.artists[i]
            mybox.set_facecolor(color_dict_soft[columns[i]])
        # change color of legend
        x_index = 0
        for xtick in bplot.get_xticklabels():
            print("xtick:%s" % xtick.get_text())
            xtick.set_color(types_color[x_index])
            x_index = x_index + 1
    # include points
    if (points):
        bplot = sns.stripplot(y=metric, x=my_column, data=emp_df, jitter=True, marker='o', alpha=0.4,   color="black",ax=ax,order=columns)

    if (annot):
        test_results = add_stat_annotation(bplot, data=emp_df, x=my_column, y=metric, box_pairs=box_pairs, order=columns, test='t-test_ind', text_format='star', loc='inside',verbose=0)

    if(metric in["EA_merged","Area_merged"]):
        ax.set_ylim(0,1.1)
    #elif(metric=="P_ratio_cLS"):
    #    ax.set_ylim(0.5, 1.1)
    #elif(metric=="cLSsim"):
    #    ax.set_ylim(0, 0.7)

    if (y_lim):
         ax.set_ylim(y_lim[0],y_lim[1])

    if (ylabel):
        label=label_dict[ylabel]
        angle=90
        if(len(label)>8):
            print("angle 90")
            angle =90
        ax.set_ylabel(label, rotation=angle, fontsize=fontsize, labelpad=12)
    else:
        print("angle 0")
        ax.set_ylabel("", rotation=0, fontsize=fontsize, labelpad=12)

    ax.set_xlabel(my_column,fontsize=21)

    if (title):
        ax.set_title('%s' % label_dict[title], fontsize=fontsize)

    if (index):
        ax.text(0.01, 0.9, index, transform=ax.transAxes, size=21, weight='bold')

    if (significance_lines):
        ax.axhline(0, color='whitesmoke', linestyle='dotted')
        ax.axhline(1.96, color='gainsboro', linestyle='dotted')
        ax.axhline(-1.96, color='gainsboro', linestyle='dotted')
        ax.axhline(2.33, color='silver', linestyle='dotted')
        ax.axhline(-2.33, color='silver', linestyle='dotted')

    return ax

def plot_regplot_fromdf_wseaborn_to_ax(emp_df,x="cLSsim",y="r_EA",reg_line=True,index=False,corr_xy=(0.05, 0.7),ax=None,**kwargs):

    fontsize = 16
    ylabel = None
    try:
        ylabel = kwargs["ylabel"]
    except:
        pass
    xlabel = None
    try:
        xlabel = kwargs["xlabel"]
    except:
        pass

    title = None
    try:
        title = kwargs["title"]
    except:
        pass

    M_df = emp_df[emp_df["sign"] != "AA"]
    A_df = emp_df[emp_df["sign"] == "AA"]

    # 1 regressions lines
    if (reg_line):
        sns.regplot(x=x, y=y, data=M_df, ax=ax, scatter=False, color=color_dict["M"], line_kws={"linewidth": 1},ci=None)
        sns.regplot(x=x, y=y, data=A_df, ax=ax, scatter=False, color=color_dict["AA"], line_kws={"linewidth": 1}, ci=None)
        #sns.regplot(x=x, y=y, data=emp_df, ax=ax, scatter=False, color='gray',                line_kws={"linewidth": 1, "linestyle": "dotted"}, ci=None)

    #scatter plots
    sns.scatterplot(x=x, y=y, hue="sign", data=emp_df, ax=ax,  palette=color_dict, marker="o", s=100, legend=False)
    corrfunc(M_df[x], M_df[y], ax, method="pearson", color=color_dict["M"], square=False,xy=(corr_xy[0], corr_xy[1]+0.1))
    corrfunc(A_df[x], A_df[y], ax, method="pearson", xy=corr_xy, color=color_dict["AA"], square=False)
    #corrfunc(emp_df[x], emp_df[y], ax, method="pearson", xy=(0.1, 0.7), square=False)
#'#008080'

    if (ylabel):
        label=label_dict[ylabel]
        angle=90
        if(len(label)>8):
            print("angle 90")
            angle =90
        ax.set_ylabel(label, rotation=angle, fontsize=fontsize, labelpad=12)
    else:
        print("angle 0")
        ax.set_ylabel("", rotation=0, fontsize=fontsize, labelpad=12)

    if (xlabel):
        ax.set_xlabel(label_dict[xlabel], rotation=0, fontsize=fontsize, labelpad=12)

    if (title):
        ax.set_title('%s' % label_dict[title], fontsize=16)


    if (index):
        ax.text(0.01, 0.9,  index, transform=ax.transAxes, size=18, weight='bold')

    return ax

def get_node_label_dictiponary(Mnet):
    sets=list(Mnet.slices[1])

    my_dict = {}

    nodes_in_set={}
    for each_set in sets:
        nodes_in_set[each_set]=get_nodes_in_set(Mnet,set_name=each_set)
        index=0
        for node in nodes_in_set[each_set]:

            my_dict[node]=each_set+"_"+str(index)
            if ("\n" in node):
                node=node.replace("\n", " ")
                my_dict[node] = each_set + "_" + str(index)
            index=index+1


    return my_dict
def plot_gml2(Mnet, K_df, name,**kwargs):

    try:
        color_by=kwargs["color_by"]
    except:
        color_by="set"

    # first obtain nodes in each set/interaction layer:
    linking_set = get_linking_set(Mnet)
    print("linking_set:%s" % linking_set)
    interactions = list(Mnet.slices[2])
    sets=list(Mnet.slices[1])
    sets.remove(linking_set)
    #K_df = get_K_df(Mnet, linking_set=linking_set)
    # nxMnet = from_pymnet_2_nx(Mnet)
    label_dict=get_node_label_dictiponary(Mnet)
    #print(label_dict)
    # read groups from Gremlin
    try:
        print("No G groups allowed\n")
        #ls_group_dict = read_Gremlin_groups(name, "")
        #print(ls_group_dict)
    except:
        print("ERROR reading groups!\n")
        ls_group_dict = {}



    y_ord = "divergent"
    try:
        y_ord=kwargs["y_ord"]
    except:
        pass

    #print(node_group)
    #xcoor_d = { linking_set : 0, sets[0]:-1, sets[1]:1 ,'flower_visitor':-1, 'leaf_miner':1,'seed_feeder':1,'butterfly':-1,'seed-feeding.insect':1, 'seed-feeding.rodent':1,'seed-feeding.bird':1,'aphid':1,'flower.visitor':-1,'Bat_f':1, 'Bat_n':-1}
    xcoor_d = {linking_set: 0, 'flower_visitor': -1, 'leaf_miner': 1, 'seed_feeder': 1, 'butterfly': -1, 'seed-feeding.insect': 1, 'seed-feeding.rodent': 1, 'seed-feeding.bird': 1, 'aphid': 1, 'flower.visitor': -1, 'Bat_f': 1, 'Bat_n': -1, 'Insect_p':1,'Insect_h':-1,'Plant':1,'Parasitoid':-1,'Floral visitor':-1,'Seed dispersal':1, 'Pollinator':-1, 'Disperser':1 , 'leaf_miner_parasitoid':-1}
    xcoor_d[linking_set]=0

    ycoor_d = {}
    for nd_set in K_df["set"].unique():
        nd_selection=K_df["set"]==nd_set
        df=K_df[nd_selection]
        if (y_ord=="O.D"):
            df["Rank"] = df["O.Degree"].rank(method="first")
            for nd in df.index:
                ycoor_d[nd]=df.loc[nd,"Rank"]/df.shape[0]

        elif (y_ord=="divergent"):
            sorted_df=df.sort_values([interactions[0], interactions[1]], ascending=[True, False])
            sorted_df.index.name="sp_name"
            sorted_df.reset_index(inplace=True)
            print("sorted_df:")
            print(sorted_df)
            ranking=sorted_df["sp_name"]
            sp_rank=pd.Series(ranking.index.values, index=ranking)
            rank=(sp_rank- sp_rank.mean())/sp_rank.shape[0]*2
            for nd in rank.index:
                ycoor_d[nd]=rank.loc[nd]

        print("nodo %s: con %s de %s, y=%s" % (nd,K_df.loc[nd,"O.Degree"], df.shape[0] ,ycoor_d[nd]))

    print(ycoor_d)

    if (color_by == "gremlin"):
        node_group = pd.Series(invert_dict_wl(ls_group_dict))
    else:
        # print(K_df)
        node_group = K_df.loc[:, "set"]
        for sp in K_df.index:
            if (K_df.loc[sp, "new_PR"] > 0):
                node_group.loc[sp] = node_group.loc[sp] + "connector"
        # print(node_group)

    filename = "../OUTPUT/Networks/Graph_%s.gdf" % (name)
    f = open(filename, "w+")
    f.write('nodedef>name VARCHAR, label VARCHAR, x DOUBLE, y DOUBLE, set VARCHAR\n')
    for node in Mnet._net.keys():
        sp_name=node[0]
        if ("\n" in sp_name):
            sp_name=sp_name.replace("\n", " ")
            print("pillada %s\n" % sp_name)
        sp_set=node[1]
        xcoor=xcoor_d[sp_set]

        #ycoor=np.random.normal(0, 1)
        ycoor=2*ycoor_d[node[0]]


        try:
            group=node_group[node[0]]
        except:
            print("sp %s NO GROUP\n" % (node[0]))
        f.write('%s, "%s", %s, %s, %s\n' % (label_dict[sp_name],sp_name,3*xcoor,ycoor,group))

    f.write('edgedef>node1 VARCHAR, node2 VARCHAR, int VARCHAR\n')
    for node in Mnet._net.keys():
        for nnode in Mnet._net[node].keys():
            #f.write("%s, %s, %s\n" % (label_dict[node[0]], label_dict[nnode[0]],node[2]))
            f.write("%s, %s, %s\n" % (label_dict[node[0]], label_dict[nnode[0]], node[2]))

    return
####### REGRESSIONS ########
def corrfunc(x, y, axi=None, method="pearson", color="black", fontsize=20,square=False, **kws):
    axi = axi or plt.gca()
    if (method=="pearson"):
        r, _ = stats.pearsonr(x, y)
    elif(method =="spearman"):
        r, _ = stats.spearmanr(x, y)

    x = 0.1
    y = 0.1
    try:
        x=kws["xy"][0]
        y=kws["xy"][1]
    except:
        pass

    try:
         color=kws["color"]
    except:
        pass

    if (square==True):
        axi.annotate(r'$r^{2} = {%.2f}$' % (pow(r,2)), xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)

    else:
        axi.annotate("r = {:.2f}".format(r),
                xy=(x, y), xycoords=axi.transAxes, color=color,fontsize=fontsize)

### HELPERS
def addlabels(x,y,label,fontsize=5,fontweight="bold",fontfamily='sans-serif'):
    for i in range(len(x)):
        plt.text(i, y[i], int(label[i]), ha = 'center',va= 'center_baseline',color="white",fontsize=fontsize,fontweight=fontweight,fontfamily=fontfamily)


#label_dict = {"Ext_Area_merged_zs_constant_NL": r'$Z_{NL}(EA)$', "Ext_Area_merged_zs_constant_Kn": r'$Z_{K}(EA)$',"EA_merged": r'$EA$', "EA_corr": r'$r_{EA}$', 'Ext_Area_corr_zs_constant_NL': r'$Z_{NL}(r_{EA})$','Ext_Area_corr_zs_constant_Kn': r'$Z_{K}(r_{EA})$', "EA_corr_RND": r'$r_{EA}(RND)$',"EA_OD": r'EA (DD)', "EA_RND": "EA (RND)", "EA_ODinv": "EA (ID)", "EA_corr_OD": r'$r_{EA}(DD)$', "EA_corr_ODinv": r'$r_{EA}(ID)$',"sigmaK": r'$\sigma_{k}$', "sigmaK_norm": r'$\sigma_{k}/<k>$', "r(k)": "r", "r_k": "r","rb_k": r'$r_{b}$',"cLSsim": "C", "PRinLhubs": r'$PR_{H}$', "CinLHubs": r'$H_{C}$', "CinLHubs_20": r'$H_{C}$', "P_ratio_cLS": r'$PR_{C}$',"PRofCinLhubs": r'$PR_{H_{C}}$', "sigmaK Zscore": r'$Z_{NL}(\sigma_{k}/<k>)$', "sigmaK_norm_Zscore": r'$Z_{NL}(\sigma_{k}/<k>)$', "r(k)_Zscore_K": r'$Z_{K}(r)$', "av_k": "<k>", "N": "N","":"","RND":"RND","OD":"DD","ODinv":"ID","Area_merged":"EA","HD":r'$\sigma_{k}/<k>$',"new_cLS_PR":r'$PR_{C}$',"LS_HD":r'$\sigma_{k}/<k>_{LS}$'}

#color_dict = {'AA': 'plum', 'MA': 'mediumaquamarine', 'MM': 'skyblue', 'H-P': '#90e4c1', 'SH-SD': '#65ab7c','H-SD': '#6fc276', 'F-N': '#008080', 'H-Pa': '#c875c4', 'SH-Pa': '#EE82EE', 'LH-SH': '#880088','P-SD': '#5a86ad', 'P-A': '#a2bffe'}

#color_dict={'AA': "#d148d3", "MA": "#615887","MM":"#519169"}
color_dict={'AA': "#663a40", "MA":  "#86b3b4","MM":"#6cbe55"}
color_dict={'AA': "#d2a088", "MA":  "#1d7583","MM":"#75be31"}
#color_dict={'AA': "#CEA0D9", "MA":  "#6DA555","MM":"#2F83BF"}
color_dict={'AA': "#d6791f", "MA": "#369094","MM":"#6cbe55"}
color_dict={'AA': "#8711ac", "MA":  "#87aa23","MM":"#4553c2", "M":"#0099cc"} #this
color_dict_soft={'AA': "#9E74A2", "MA":  "#87aa23","MM":"#828BD1"}

    #renaming dict for all metrics ever used
label_dict={"Ext_Area_merged_zs_NL2": r'$Z_{NL2}(EA)$', "Ext_Area_merged_zs_K2": r'$Z_{K2}(EA)$',"Ext_Area_merged_zs_NL": r'$Z_{NL}(EA)$', "Ext_Area_merged_zs_K": r'$Z_{K}(EA)$', "EA_merged": r'$EA$', "R_merged":r'$R_{M}$', "EA_corr":r'$r_{EA}$' , 'Ext_Area_corr_zs_NL': r'$Z_{NL}(r_{EA})$', 'Ext_Area_corr_zs_K': r'$Z_{K}(r_{EA})$',"sigmaK_norm": r'$\sigma_{k}/<k>$' ,"sigmaK_norm_Zscore": r'$Z_{NL}(\sigma_{k}/<k>)$ ',"r(k)":"r", "r(k)_Zscore_NL": r'$Z_{NL}(r)$', "r(k)_Zscore_K": r'$Z_{K}(r)$',"k_sigma_norm":r'$\sigma_{k}/<k>$' ,"K_sig. Zscore_norm":r'$Z(\sigma_{k}/<k>)$ NL',"Z_r_NL":r'$Z(r)$ NL',"Z_r_K":r'$Z(r)$ K',"epsi":r'$\epsilon$', "epsiZscore":r'$Z(\epsilon)$', "r":"r", "cLSsim":"C","CinLHubs":r'$H_{c}$',"P_ratio_cLS":r'$PR_{C}$','Ext_Area_corr_zs_NL2':r'$Z_{NL2}(r_{EA})$','Ext_Area_corr_zs_K2':r'$Z_{K2}(r_{EA})$', 'Ext_Area_corr_zs_constant_KnonLS_free':r'$Z_{F}(r_{EA})$', "Ext_Area_merged_zs_constant_KnonLS_free": r'$Z_{F}(EA)$', "HinC_10pc":r'$C_{H_{10}}$',"HinC_5":r'$C_{H_{5}}$',"new_PR":r'$PR_{LS}$',"new_cLS_PR": r'$PR_{C}$','rb(k)':r'$r_{B}$','rb_k':r'$r_{B}$','r_clean':r'$\rho$', "emp_LS_HD":r'$\sigma_{k_{LS}}/<k_{LS}>$',"emp_nonLS_HD":r'$\sigma_{k_{nLS}}/<k_{nLS}>$',"LS_HD_ZS":r'$\sigma_{k_{LS}}/<k_{LS}>$', "nonLS_HD_ZS":r'$\sigma_{k_{nLS}}/<k_{nLS}>$', "P_HD":r'$\sigma_{k_{P}}/<k_{P}>$', "nonP_HD":r'$\sigma_{k_{nP}}/<k_{nP}>$', "emp_P_HD":r'$\sigma_{k_{P}}/<k_{P}>$',"emp_nonP_HD":r'$\sigma_{k_{nP}}/<k_{nP}>$', 'P_HD_ZS':r'$\sigma_{k_{P}}/<k_{P}>$', 'nonP_HD_ZS':r'$\sigma_{k_{nP}}/<k_{nP}>$',"":"","RND":"RND","OD":"DD","ODinv":"ID", "HD":r'$\sigma_{k}/<k>$', "LS_HD":r'$\sigma_{k_{LS}}/<k_{LS}>$', "CinLHubs_20": r'$H_{C}$',"r_k": "r",'Zs_NL_HD':r'$Z_{NL}(\sigma_{k}/<k>)$','Zs_K_r_k': r'$Z_{K}(r)$', 'Zs_K_rb_k': r'$Z_{K}(r_b)$','Zs_NL_r_EA':r'$Z_{NL}(r_{EA})$', 'Zs_K_r_EA':r'$Z_{K}(r_{EA})$','Zs_NL2_r_EA':r'$Z_{NL2}(r_{EA})$','Zs_K2_r_EA':r'$Z_{K2}(r_{EA})$','Zs_NL2_Area_merged':r'$Z_{NL2}(EA)$','Zs_K2_Area_merged': r'$Z_{K2}(EA)$','Zs_NL_Area_merged':r'$Z_{NL}(EA)$','Zs_K_Area_merged':r'$Z_{K}(EA)$', "Area_merged":r'$EA_{M}$','lsB_HD':r'$(\sigma_{k}/<k>)_{LS_B}$','lsA_HD':r'$(\sigma_{k}/<k>)_{LS_A}$','LS_HD':r'$(\sigma_{k}/<k>)_{LS}$'}