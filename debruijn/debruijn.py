#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "BAILLIF Marine"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Marine BAILLIF"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    with open(fastq_file,"r") as fastq : 
        for ligne in fastq :
            yield next(fastq).strip()
            next(fastq)
            next(fastq)
    pass



def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    for i in range(len(read) - (kmer_size-1)):
        yield (read[i:i+kmer_size])
    pass



def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    list_seq = list(read_fastq(fastq_file))
    list_read = []
    for i in range(len(list_seq)):
        sequence = list_seq[i]
        list_kmer = list(cut_kmer(sequence, kmer_size))
        list_read.append(list_kmer)
        #list_read.append(list_kmer)
    
    #print(list_read)
    dico_kmer = {}
    for i in range(len(list_read)):
        for kmer in list_read[i]:
            if kmer in dico_kmer:
                dico_kmer[kmer] += 1
            else:
                dico_kmer[kmer] = 1
    return(dico_kmer)
    pass


def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    G=nx.DiGraph()
    for kmer in kmer_dict:
        kmer_all= kmer
        
        kmer_start = kmer_all[:-1]
        kmer_end = kmer_all[1:]
        
        if kmer_start not in G:
            G.add_node(kmer_start)
        if kmer_end not in G:
            G.add_node(kmer_end)
        
        G.add_edge(kmer_start,kmer_end,weight=kmer_dict[kmer_all])
    return G
    pass



def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    #list_no_pred = []
    #for node in graph : 
    #    print(node)
    #    pred_nod = list(nx.predecessor(graph, node))
    #    if len(pred_nod) == 0 :
    #        list_no_pred.append(pred_nod)
    #return list_no_pred

    nodes = graph.nodes()
    list_no_pred = []

    for node in nodes:
        pred_node=list(graph.predecessors(node))
        if len(pred_node)==0:
            list_no_pred.append(node)
    return list_no_pred
    pass

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """

    nodes = graph.nodes()
    list_no_succ = []

    for node in nodes:
        success_node = list(graph.successors(node))
        if len(success_node)==0:
            list_no_succ.append(node)
    return list_no_succ
    pass

def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    list_contig = []
    for node_start in starting_nodes:
        for node_end in ending_nodes:
            if nx.has_path(graph, node_start, node_end) == True : 
                print("list_seq")
                liste_seq = list(nx.all_simple_paths(graph, node_start,node_end))
            #print("contig")
                for liste in liste_seq : 
                    seq = liste[0]
                    for mini_seq in liste[1:] :
                        seq = seq + mini_seq[-1]
                    #print((seq,len(seq)))
                    list_contig.append((seq,len(seq)))
                
            #if len(liste_seq) != 0:
                #seq = liste_seq[0]
                #print(f"seq {seq}")
                #for mini_seq in liste_seq :
                    #print("mini_seq") 
                    #print(mini_seq)
                    #print(mini_seq[-1])
                   #seq =  seq + mini_seq[-1]
                #list_contig.append((seq,len(seq)))
    print(list_contig)
    return list_contig
    pass

def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    pass


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    diko_kmer = build_kmer_dict("eva71_two_reads.fq",22)
    #print(diko_kmer)

    graph_kmer = build_graph(diko_kmer)

    list_predecessor = get_starting_nodes(graph_kmer)
    #print(list_predecessor)

    #sucessor
    list_sucessor = get_sink_nodes(graph_kmer)
    #print(list_sucessor)

    # contig 
    contig = get_contigs(graph_kmer,list_predecessor,list_sucessor)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == '__main__': # pragma: no cover
    main()

# read fastq
#mon_generateur = list(read_fastq("eva71_two_reads.fq"))
#for i in mon_generateur : 
    #print(i)

# cut kmer 
#generateur2 = cut_kmer(mon_generateur[0],3)
#for i in generateur2 : 
    #print(i)