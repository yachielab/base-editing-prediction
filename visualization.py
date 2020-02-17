import os 
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ["Helvetica","Arial","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['figure.figsize']  = [3,3]
matplotlib.rcParams['font.size']         = 15
matplotlib.rcParams["axes.labelpad"]     = 10.0 
matplotlib.rcParams["axes.labelcolor"]   = "#000000"
matplotlib.rcParams["axes.linewidth"]    = 0.8 
matplotlib.rcParams["xtick.major.width"] = 0.8
matplotlib.rcParams["xtick.major.pad"]   = 4.0
matplotlib.rcParams["ytick.major.width"] = 0.8
matplotlib.rcParams["xtick.major.pad"]   = 4.0
cmap1 = plt.cm.tab10
cmap2 = plt.cm.Set3  
colors1 = [cmap1(i) for i in range(0,10)]
colors2 = [cmap2(i) for i in range(0,12)] 

color_dict  = {"G":"#f2f059", "C":"#74b2d7", "A":"#79E5B7", "T":"#ff776c", "N":"#FFFFFF", "-":"#FFFFFF"}

def main(ref,editing_spec,s,e,name):
    spec = [] 
    fig  = plt.figure(figsize=(4,2))
    
    #Visualization for the editing spectrum 
    ax   = fig.add_axes([0.1,0.2,0.8,0.8])
    for i in range(len(ref)):
        bottom = 0
        for c in ["A","T","G","C"]:
            value = editing_spec[i][c]
            if value != "N.A.":
                ax.bar([i], [value], width=0.9, align="center", bottom=bottom, facecolor=color_dict[c])
                bottom += editing_spec[i][c]
        spec.append(bottom) 

    ax.spines["top"].set_visible(False) 
    ax.spines["right"].set_visible(False)
    ax.set_xlim(-0.5,len(ref)-0.5)
    if max(spec) > 0.5:
        ax.set_ylim(0,1.0) 
    elif max(spec) > 0.3: 
        ax.set_ylim(0,0.5) 
    else: 
        ax.set_ylim(0,0.3) 
    ax.set_xticks([]) 
    ax.set_ylabel("Substitution frequency")
    ax.set_title(name,fontsize=15) 
    
    #Visualization for reference sequence
    ax  = fig.add_axes([0.1,0.13,0.8,0.07])
    for i in range(len(ref)):
        ax.bar([i], [1], bottom=0, width=1.0, lw=0.2, edgecolor="#BBBBBB", facecolor=color_dict[ref[i]])
    ax.spines["top"].set_visible(False) 
    ax.spines["bottom"].set_visible(False) 
    ax.spines["left"].set_visible(False) 
    ax.spines["right"].set_visible(False)
    ax.set_xlim(-0.5,len(ref)-0.5)
    ax.set_ylim(0,1.0) 
    ax.set_yticks([])
    ax.set_xticks([]) 

    #Visualization for PAM and target sequence
    ax  = fig.add_axes([0.1,0.08,0.8,0.05])
    for i in range(len(ref)):
        if -20 <= i + s <= -1: 
            ax.bar([i], [1], bottom=0, width=1.0, lw=0.0, edgecolor="#BBBBBB", facecolor="#808080")
        elif 0 <= i + s <= 2:
            ax.bar([i], [1], bottom=0, width=1.0, lw=0.0, edgecolor="#BBBBBB", facecolor="#FF2222")
        else:
            ax.bar([i], [1], bottom=0, width=1.0, lw=0.0, edgecolor="#DDDDDD", facecolor="#EEEEEE")
    pos = 0 - s
    if e < 0:
        ticks = [pos-20,pos-10,pos-1] 
        ax.set_xticklabels(["-20","-10","-1"]) 
    elif e < 2:
        ticks = [pos-20,pos-10,pos] 
        ax.set_xticklabels(["-20","-10","0"]) 
    else: 
        ticks = [pos-20,pos-10,pos,pos+2]
        ax.set_xticklabels(["-20","-10","0","+2"]) 

    ax.spines["top"].set_visible(False) 
    ax.spines["bottom"].set_visible(False) 
    ax.spines["left"].set_visible(False) 
    ax.spines["right"].set_visible(False) 
    ax.set_xticks(ticks) 
    ax.set_xlim(-0.5,len(ref)-0.5)
    ax.set_ylim(0,1.0) 
    ax.set_yticks([])
    
    #Legends for PAM and target sequence
    if e > -1:
        ax_pam    = fig.add_axes([-0.05,-0.2,0.03,0.06]) 
        ax_target = fig.add_axes([0.1,-0.2,0.03,0.06])
        for ax, char, color in zip((ax_pam, ax_target),("PAM","Target"),("#FF0000","#808080")):  
            ax.set_xticks([]) 
            ax.set_yticks([]) 
            ax.spines["top"].set_linewidth(0.0)
            ax.spines["bottom"].set_linewidth(0.0)
            ax.spines["right"].set_linewidth(0.0)
            ax.spines["left"].set_linewidth(0.0)
            ax.patch.set_facecolor(color)
            ax.text(1.3,0.4,char,rotation=0,va="center",ha="left",fontsize=13) 
    else: 
        ax = fig.add_axes([0.1,-0.2,0.03,0.06])
        ax.set_xticks([]) 
        ax.set_yticks([]) 
        ax.spines["top"].set_linewidth(0.0)
        ax.spines["bottom"].set_linewidth(0.0)
        ax.spines["right"].set_linewidth(0.0)
        ax.spines["left"].set_linewidth(0.0)
        ax.patch.set_facecolor("#808080")
        ax.text(1.3,0.4,"Target",rotation=0,va="center",ha="left",fontsize=13) 

    #Legends for nuclotides
    ax_A = fig.add_axes([0.1 + 0.4,-0.2,0.03,0.06]) 
    ax_T = fig.add_axes([0.1 + 0.5,-0.2,0.03,0.06]) 
    ax_G = fig.add_axes([0.1 + 0.6,-0.2,0.03,0.06]) 
    ax_C = fig.add_axes([0.1 + 0.7,-0.2,0.03,0.06]) 
    for ax, char in zip((ax_A, ax_T, ax_G, ax_C),("A","T","G","C")):  
        ax.set_xticks([]) 
        ax.set_yticks([]) 
        ax.spines["top"].set_linewidth(0.0)
        ax.spines["bottom"].set_linewidth(0.0)
        ax.spines["right"].set_linewidth(0.0)
        ax.spines["left"].set_linewidth(0.0)
        ax.patch.set_facecolor(color_dict[char])
        ax.text(1.3,0.4,char,rotation=0,va="center",ha="left",fontsize=13)

    fig.savefig("editing_spectrum.pdf",bbox_inches="tight") 
