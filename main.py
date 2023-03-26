import networkx as nx
import networkx.algorithms.community as nx_comm
import copy
import csv

# Data Representation
# 1. CS=[[1,4,2],[2,3,5]]
# 2. SubG=netowrkx object

global_graph=nx.Graph()

# ---------------------------------------------------------------
# Algorthim 3:  Identifying the types of incremental elements

def isSubset(V, C):
    for i in range(len(V)):
        if(V[i] not in C):
            return False
    return True

def incremental_elements(CSt, SubG):
    number_of_communities = len(CSt)
    Vold=[]
    RelComLab=set()
    for j in range(number_of_communities):
        Voldj=list(set(list(SubG.nodes))&set(CSt[j]))
        if(set(Voldj)==set(list(SubG.nodes))):
            type=2
            RelComLab.add(j)
            return type, RelComLab
        if(len(Voldj)!=0):
            Vold.append(Voldj)
            RelComLab.add(j)
    n=len(Vold)
    if n==0:
        type=1
    if n==1:
        type=3
    if n>1:
        type=4
    return type, RelComLab
# ---------------------------------------------------------------


# ---------------------------------------------------------------
# Algorithm 2: Update communities for the MultiCont type.

def calculate_subordinating_stength_subg(v,graph):
    neighbors=list(global_graph.neighbors(v))
    members=list(graph.nodes())
    denom=0
    numer=0
    for node in neighbors:
        denom+=global_graph[node][v]['weight']
        if node in members:
            numer+=global_graph[node][v]['weight']
    return numer/denom

def calculate_subordinating_stength_community(v,community):
    neighbors=list(global_graph.neighbors(v))
    denom=0
    numer=0
    for node in neighbors:
        denom+=global_graph[node][v]['weight']
        if node in community:
            numer+=global_graph[node][v]['weight']
    return numer/denom

def update_communities_multicont(CS,SubG):
    k=len(CS)
    V=[[] for _ in range(k)]
    for i in range(k):
        V[i]=list(set(list(SubG.nodes))&set(CS[i]))
    V_set=set()
    for i in range(k):
        for node in V[i]:
            V_set.add(node)
    for node in V_set:
        max_s_strength=0
        index=0
        for i in range(k):
            if node in CS[i]:
                s_strength=calculate_subordinating_stength_community(node,CS[i])
                if s_strength>max_s_strength:
                    max_s_strength=s_strength
                    index=i
        if max_s_strength>calculate_subordinating_stength_subg(node,SubG):
            SubG.remove_node(node)
            CS[index].append(node)
    updateCS=CS
    updateCS.append(list(SubG.nodes))
    return updateCS
# ---------------------------------------------------------------


# ---------------------------------------------------------------
# Algorithm 1 Update communities for the Mixed type.

def update_communities_mixed(Vold,Vnew,Ck,SubG):
    updateCS=[]
    for node in Vold:
        if calculate_subordinating_stength_community(node,Ck)>calculate_subordinating_stength_subg(node,SubG):
            SubG.remove_node(node)
            if node not in Ck:
                Ck.append(node)
    for node in Vnew:
        if calculate_subordinating_stength_community(node,Ck)>calculate_subordinating_stength_subg(node,SubG):
            SubG.remove_node(node)
            if node not in Ck:
                Ck.append(node)
    updateCS.append(Ck)
    updateCS.append(list(SubG.nodes))
    return updateCS
    
# ---------------------------------------------------------------


# ---------------------------------------------------------------
# Algorithm 4 Detecting communities in the evolving network incrementally
def update_edge_weight(graph,alpha,theta):
    nodes=list(global_graph.nodes)
    for u in nodes:
        for v in nodes:
            if not global_graph.has_edge(u,v):
                continue
            temp= graph[u][v]["weight"] if graph.has_edge(u,v) else 0
            temp+=alpha*global_graph[u][v]["weight"]
            global_graph[u][v]["weight"]=temp if temp>theta else 0

def modularity(Ckt):
    subg=global_graph.subgraph(Ckt).copy()
    orgmod=nx_comm.modularity(subg,[set(Ckt)])
    edges=list(global_graph.edges.data("weight", default=1))
    edges.sort(key = lambda x: x[2])
    for edge in edges:
        demo_graph=copy.deepcopy(global_graph)
        demo_graph.remove_edge(edge[0],edge[1])
        mod=nx_comm.modularity(demo_graph,Ckt)
        if (mod-orgmod) > 0:
            subg.remove_edge(edge[0],edge[1])
            global_graph.remove_edge(edge[0],edge[1])
    
    return [list(subg.subgraph(c).copy().nodes()) for c in nx.connected_components(subg)]

def get_communities_from_incremental_elements(CSt,delta_G,alpha,theta):
    Set_subG = [delta_G.subgraph(c).copy() for c in nx.connected_components(delta_G)]
    uplabel=set()
    CSt1=[]
    update_edge_weight(delta_G,alpha,theta)
    for SubG in Set_subG:
        type,RelComLab=incremental_elements(CSt,SubG)
        print(type)
        uplabel=uplabel&RelComLab
        if type==1:
            CSt1.append(list(SubG.nodes))
        
        if type==2:
            for rel in RelComLab:
                newcommunties=modularity(CSt[rel])
                CSt[RelComLab]=newcommunties[0]
                for i in range(1,len(newcommunties)):
                    CSt.append(newcommunties[i])
                    uplabel.add(CSt.index(newcommunties[i]))

        if type==3:
            for rel in RelComLab:
                Vold=list(set(list(SubG.nodes))|set(CSt[rel]))
                Vall=list(SubG.nodes)
                Vnew = list(set(Vall)-set(Vold))
                newcommunties=update_communities_mixed(Vold,Vnew,CSt[rel],SubG)
                CSt[rel]=newcommunties[0]
                for i in range(1,len(newcommunties)):
                    CSt.append(newcommunties[i])
                    uplabel.add(CSt.index(newcommunties[i]))
        if type==4:
            CSt=update_communities_multicont(CSt,SubG)
    for i in range(len(CSt)):
        if i not in uplabel:
            newcommunties=modularity(CSt[i])
            for j in range(0,len(newcommunties)):
                CSt1.append(newcommunties[j])
        else:
            CSt1.append(CSt[i])
    return CSt1

# ---------------------------------------------------------------

with open("initialgraph.csv", mode ='r')as file:
    csvFile = csv.DictReader(file)
    for line in csvFile:
        global_graph.add_edge(int(line['source']),int(line['target']),weight=int(line['weight']))
CS=[[0,1,2,3],[4,5,6],[7,8,9],[10,11,12,13,14,15]]

subG=nx.Graph()
with open("completed_contained.csv", mode ='r')as file:
    csvFile = csv.DictReader(file)
    for line in csvFile:
        subG.add_edge(int(line['source']),int(line['target']),weight=int(line['weight']))
        if not global_graph.has_edge(int(line['source']),int(line['target'])):
            global_graph.add_edge(int(line['source']),int(line['target']),weight=int(line['weight']))
        else:
            global_graph[int(line['source'])][int(line['target'])]['weight']+=int(line['weight'])

Cst=get_communities_from_incremental_elements(CS,subG,0.5,3.5)
print(Cst)