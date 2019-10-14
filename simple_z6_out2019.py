#python simple-z6_ana_set.2019.py K noise nsamples N T
# -*- coding: utf-8 -*-
import numpy as np
from numpy import ndarray
#import matplotlib.pyplot as m
import math as math
import networkx as nx
import cmath as cm
import sys

class node:
    def __init__(self, zcomplex, ang_freq, neighbors):
        self.z = zcomplex
#        self.b = b
        self.w = ang_freq
        self.n = neighbors

    def f(self):
        return (a*abs(self.z)**4+b*abs(self.z)**2+c+self.w*1j)*self.z

 #   def db(self):
 #       return (self.bo * (1 - abs(self.z) /self.zo) - self.b) / tau

    def fut(self):
        self.z+=self.f()*dt
#        self.b+=self.db()*dt

    def mov(self):
        soma=0.0
        for i,w in enumerate(self.n):
            soma+=nod[w].z
        soma=beta*soma*dt
        self.z+=noise*(np.random.normal(0,1.)+np.random.normal(0,1.)*1.j)*math.sqrt(dt)+soma/(N-1.)

#Function to plot covariance and average histogram (set printhisto=1)

def histograma(abscov,absav,absicto,count):
    
    axs[0][count].hist(abscov,bins=40, range=[0,1],label=str(p), stacked=True, fill=False)
    axs[0][count].set_title(str(p))
    axs[1][count].hist(absav,bins=40, range=[0,1], stacked=True, fill=False)
   # axs[1][count].set_title(str(p))
    axs[2][count].hist(absicto,bins=40, range=[0,1], stacked=True, fill=False)
   # axs[2][count].set_title(str(p))
    if(count==3):
        title="N=" + str(N) + " T="+ str(T)  + "Correlation (up)  Averages (middle)  Ictogenicity (down)"
        m.suptitle(title)
        print axs
    return axs

# Function to plot nodes neighbor list (set printnet=1)

def netprint(node_neigh_list):
    for i,w in enumerate(node_neigh_list):
        rede.write("{:d} - ".format(i))
        for j in w:
            rede.write("{:d} ".format(j))
        rede.write("\n")
               
# Here starts the main program       
entrada = sys.argv[0:6]
K=int(entrada[1])
global a,b,beta,c,dt,noise#,N
a=-1
b=2.
c=-0.9
icto=math.sqrt((-b+math.sqrt(b**2-4*a*c))/(2*a))
beta=1.4
dt=0.01
noise=float(entrada[2])
print K,noise
nsamples = int(entrada[3])
N=int(entrada[4])#100
T=int(entrada[5])#2500
Nij=N*(N-1)/2
printhisto=0 #set to 1 to print the histogram
printnet=1 #set to 1 to print the nodes neighbor list
#fig, axs = m.subplots(nrows=3,ncols=2)

count=0
plist=[0.01,0.1,0.5,0.6,0.7,0.8,0.9,1.0]
meanavicto = ndarray((len(plist),),float)
desviocov = ndarray((len(plist),),float)
meancov = ndarray((len(plist),),float)

print(plist.index(0.01))

for p in plist:
    meanavicto[plist.index(p)] = 0.0
    desviocov[plist.index(p)] = 0.0
    meancov[plist.index(p)] = 0.0
    
saida2 = open("avsample_z6_N" + str(N) + "_T" + str(T) + "_k" + str(K) + "_d" + str(noise) + ".out","w")
saida2.write("#p, cor. media, desvio cor. media, ictogenicity \n")
for sample in range(1,nsamples+1):
    saida=open("z6_N" + str(N) + "_T" + str(T) + "_k" + str(K) + "_d" + str(noise) + "_sample_"+str(sample)  +".out","w")
    saida.write("#p, cor. media, desvio cor. media, ictogenicity \n")
    for p in plist:
        #for p in [0.9, 1.0]:
        print 'p=',p
        # pp=int(math.floor(p*(N-1)))
        rede=open("z6-net-N" + str(N) + "k" + str(K) + "_d" + str(noise) + "_sample_"+str(sample)  + ".out","w")
        rede.write("p= {:f}\n".format(p))
        #    bo=5
        #    bth=2*math.sqrt(a*c)
        #    tau=1000
        #    zo=-0.5*bth/a
        #Constructing connection network
        g=nx.Graph()
        g.add_node(0)
        g.add_node(1)
        while not nx.is_connected(g):
            g = nx.watts_strogatz_graph(N,K,p) # número de nós, grau médio, probabilidade de reconexão
        #        g = nx.barabasi_albert_graph(N,m) # número de nós, número de vizinhos de cada nó novo
        
        node_neigh_list=[]
        for i in nx.nodes(g):
            list_neigh=[]
            for j in nx.all_neighbors(g, i):
                list_neigh.append(j)
            node_neigh_list.append(list_neigh)
                
                
        if(printnet==1):
            netprint(node_neigh_list)
                    
        #initialize N nodes
        nod=list(node(np.random.normal(0,.1)+np.random.normal(0,.1)*1.j, np.random.uniform(-0.2,0.2), node_neigh_list[i]) for i in range(N))
                    
        #Node evolution
        t=0
        x=[]
        d=[]
        #b=[]
        k=0
#        meanavicto =0.0
        avicto=0.0
#        desviocov = 0.0
#        meancov = 0.0
        cov=list(0 for i in range(Nij))
        avphi=list(0 for i in range(N))
        while(t<T):
            t+=dt
            # if(t%10<dt):
            #print "t=",int(t)
            map(lambda i:i.fut(), nod)
            map(lambda i:i.mov(), nod)
        
            #    x.append(t)
            k=0
            for i in range(N):
                avphi[i]+=cm.exp(1j*(cm.phase(nod[i].z)))
                if(abs(nod[i].z)>icto):
                    avicto+=1.
                for j in range(i+1,N):
                    cov[k]+=dt*cm.exp(1j*(cm.phase(nod[i].z)-cm.phase(nod[j].z)))
                    #                cov[k]+=-dt*cm.exp(1j*(cm.phase(nod[i].z)))*dt*cm.exp(1j*(-cm.phase(nod[j].z)))/T
                    k+=1
        #d.append(map(lambda i:cm.phase(i.z), nod))
        #   b.append(map(lambda i:i.b, nod))
                
                
        #x=list(i for i in range(Nij))
        abscov=map(lambda i:abs(i)/T, cov)
        absav=map(lambda i:abs(i)*dt/T, avphi)
        avicto=avicto*dt/T/N
        meanavicto[plist.index(p)]+=avicto
        absicto=avicto
        if(printhisto==1):
            histograma(abscov,absav,absicto,count)
            count+=1
        
        avabscov=sum(abscov)/Nij
        meancov[plist.index(p)]+=avabscov
        desviocov[plist.index(p)]+=avabscov*avabscov
        desvio_correlacao=sum(map(lambda i:(i-avabscov)**2, abscov))/Nij
        saida.write('{:f} {:f} {:f} {:f}\n'.format(p,avabscov,desvio_correlacao,avicto))
        
meanavicto = meanavicto/nsamples
meancov = meancov/nsamples
for i in range(len(plist)):
    desviocov[i] = math.sqrt(desviocov[i]/nsamples-meancov[i]*meancov[i])
for p in plist:
    saida2.write('{:f} {:f} {:f} {:f}\n'.format(p,meancov[plist.index(p)],desviocov[plist.index(p)],meanavicto[plist.index(p)]))
#Graph
#    m.subplot(2,1,1)
#m.plot(x,b)
#    m.subplot(2,1,2)
#m.plot(x,abscov)
#name = "dados_k=" + str(K) + "_d=" + str(noise) +  ".dat"
#arquivo_saida=open(name,'w')
#for i,w in enumerate(axs) :
#    f.write('%d %f',(i,w))
#fig.savefig(name,dpi=300)
#m.show()
