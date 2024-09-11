import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
import matplotlib

#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from itertools import chain, combinations



# # Add color to edges
def add_color(set_DDIs):
    ColorLegend = {}

    col = colors.cnames
    color = list(col.values())

    # --------------remove colors-------------
    index = color.index('#FFEBCD')
    color.pop(index)
    index = color.index('#FF7F50')
    color.pop(index)
    index = color.index('#DEB887')
    color.pop(index)
    index = color.index('#FFF8DC')
    color.pop(index)
    index = color.index('#DC143C')
    color.pop(index)
    # --------------remove colors-------------

    effect_impact = list(pd.value_counts(set_DDIs.effect_impact).index)

    for i in range(len(effect_impact)):
        set_DDIs.loc[set_DDIs.effect_impact == effect_impact[i], "edge_color"] = color[i + 8]
        ColorLegend[effect_impact[i]] = color[i + 8]

    set_DDIs.reset_index(inplace=True)
    set_DDIs = set_DDIs.drop(columns=['index'])
    return color, ColorLegend, set_DDIs


# # Visualize graph
def add_di_edge_to_graph(set_ddi, g, dict_control):
    # Add edges and edge attributes
    for i, elrow in set_ddi.iterrows():
        g.add_edge(elrow[0], elrow[1], attr_dict=elrow[3], create_using=nx.MultiDiGraph())
        if (elrow[0], elrow[1]) in dict_control:
            dict_control[(elrow[0], elrow[1])] = dict_control[(elrow[0], elrow[1])] + 1
        else:
            dict_control[(elrow[0], elrow[1])] = 1
    return g, dict_control


def edge_color(g):
    # Define data structure (list) of edge colors for plotting
    edge_colors = [e[2]['attr_dict'] for e in g.edges(data=True)]
    return edge_colors


def get_multiple_edge(g, dict_control):
    multiple_edge = []
    count = 0
    pair = ('', '')
    for e in g.edges:
        if dict_control[(e[0], e[1])] == 1:
            multiple_edge.append(0.0)
        elif pair != (e[0], e[1]):
            count = dict_control[(e[0], e[1])]
            for i in range(count):
                multiple_edge.append(0.06 * (i + 1))
            pair = (e[0], e[1])

    return multiple_edge


def get_node_size(g):
    list_node = list(g.nodes)
    d2_size = []
    for n2 in list_node:
        d2_size.append(150 + (g.in_degree(n2) * 50))
    return d2_size


# # Color nodes
def get_node_color(g, label_dsd, color_mark, color_basic):
    list_node = list(g.nodes)
    print('check:', list_node)
    n_color = [color_basic] * g.number_of_nodes()
    for d in label_dsd:
        index = list_node.index(d)
        n_color[index] = color_mark
    return n_color


# # Position of nodes in bipartite graph
def get_position_nodes(g, set_dsd_label):
    bottom_nodes, top_nodes = bipartite.sets(g)
    print('top_nodes:', top_nodes, 'set_dsd_label:', set_dsd_label)
    pos = dict()
    if top_nodes != set_dsd_label:
        temporal_right = top_nodes
        top_nodes = bottom_nodes
        bottom_nodes = temporal_right
    if len(bottom_nodes) > len(top_nodes):
        gap = len(bottom_nodes) / len(top_nodes)
        pos.update((n, (0.0, i)) for i, n in enumerate(bottom_nodes))
        pos.update((n, (0.005, (i * gap))) for i, n in enumerate(top_nodes))
    else:
        gap = len(top_nodes) / len(bottom_nodes)
        pos.update((n, (0.0, i * gap)) for i, n in enumerate(bottom_nodes))
        pos.update((n, (0.005, (i))) for i, n in enumerate(top_nodes))
    return pos, bottom_nodes, top_nodes


def plot_bipartite_graph(g, n_color, ColorLegend, d2_size, set_dsd_label, edge_colors, multiple_edge, plot_name, dim):
    plt.figure(figsize=dim)
    ax = plt.subplot(1, 1, 1)

    pos, bottom_nodes, top_nodes = get_position_nodes(g, set_dsd_label)
    nx.draw_networkx_nodes(g, pos=pos, node_color=n_color, node_size=d2_size)  # 'skyblue'  , label=True
    print(pos)

    i = 0
    for edge in g.edges:
        if edge[0] in list(set_dsd_label):
            style = "<|-"
        else:
            style = "-|>"
        # print(edge, style)
        # ax.annotate("", xy=pos[edge[0]], xytext=pos[edge[1]], arrowprops=dict(arrowstyle="<|-", color='moccasin'))

        ax.annotate("", xy=pos[edge[0]], xytext=pos[edge[1]],
                    arrowprops=dict(arrowstyle="<|-, head_length=0.9, head_width=0.6",
                                    color=edge_colors[i], lw=1.5,
                                    connectionstyle="arc3,rad=" + str(multiple_edge[i]), ), )
        i += 1

    nx.draw_networkx_labels(g, pos, font_size=10, alpha=1.0)  # ,font_color='r'  , font_weight="bold"

    # ax.margins(0.2, 0.2)
    ax.set_xmargin(0.35)
    plt.axis('off')
    plt.title('Drug-Drug Interactions', fontsize=12, ha='center')

    for label in ColorLegend:
        ax.plot([0], [0], color=ColorLegend[label], label=label)  # , marker='o'
    plt.legend(bbox_to_anchor=(0.5, 0), loc='upper center', fontsize=10, ncol=2).get_frame().set_alpha(0.0)
    plt.tight_layout()

    plt.savefig(plot_name + '.pdf')  # _optimal
    plt.show()
    return bottom_nodes, top_nodes


def run_plot_graph(g, n_color, ColorLegend, d2_size, edge_colors, multiple_edge, use_case, dim):
    plt.figure(figsize=dim)#figsize=(16, 8)
    ax = plt.subplot(1,1,1)
    
    pos = nx.circular_layout(g, center=(1.06,0.0), scale=1)# scale=1
    nx.draw_networkx_nodes(g, pos=pos, node_color=n_color, node_size=d2_size)

    i=0
    for edge in g.edges:
        ax.annotate("", xy=pos[edge[0]], xytext=pos[edge[1]],
                    arrowprops=dict(arrowstyle="<|-, head_length=0.9, head_width=0.6",
                                    color=edge_colors[i], lw=1.5,
                                    connectionstyle="arc3,rad="+str(multiple_edge[i]), ), )
        i+=1

    nx.draw_networkx_labels(g, pos, font_weight="bold", alpha=1.0, font_size=10) #, font_size=10
    ax.set_xmargin(0.06)
    
    plt.axis('off')
    
    for label in ColorLegend:
        ax.plot([0],[0], color=ColorLegend[label], label=label)#, marker='o'
    plt.legend(loc='upper left', fontsize=10, ncol=3, bbox_to_anchor=(0.0, 0.2)).get_frame().set_alpha(0.0)#upper center bbox_to_anchor=(0.5, 1.08)
    plt.tight_layout()
    
    use_case = use_case.replace('\n','_')
    plt.savefig(use_case+'.pdf')
    #plt.show()
    

def preprocess(union_plot, set_dsd_label):
    #union = combine_col(union)
    color, ColorLegend, union = add_color(union_plot)

    g = nx.MultiDiGraph()  # DiGraph
    dict_control = {}

    g, dict_control = add_di_edge_to_graph(union, g, dict_control)
    edge_colors = edge_color(g)
    multiple_edge = get_multiple_edge(g, dict_control)

    n_color = get_node_color(g, list(set_dsd_label), 'red', 'skyblue')
    d2_size = get_node_size(g)
    
    return g, n_color, ColorLegend, d2_size, set_dsd_label, edge_colors, multiple_edge
    
    
    
def plot_bipartite(union_plot, set_dsd_label, plot_name):
    g, n_color, ColorLegend, d2_size, set_dsd_label, edge_colors, multiple_edge = preprocess(union_plot, set_dsd_label)
    bottom_nodes, top_nodes = plot_bipartite_graph(g, n_color, ColorLegend, d2_size, set_dsd_label, edge_colors,
                                                   multiple_edge, plot_name, (6, 5))  # (8,6)
    
    
    
def plot_graph(union_plot, set_dsd_label, plot_name):
    g, n_color, ColorLegend, d2_size, set_dsd_label, edge_colors, multiple_edge = preprocess(union_plot, set_dsd_label)
    run_plot_graph(g, n_color, ColorLegend, d2_size, edge_colors, multiple_edge, plot_name, (8, 5)) #(12, 6)
    
    
