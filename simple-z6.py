import numpy as np
import matplotlib.pyplot as m
import math as math
import networkx as nx
import cmath as cm

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

def histograma(abscov,absav,count):
    axs[0][count].hist(abscov,bins=40, range=[0,1],label=str(p), stacked=True, fill=False)
    axs[0][count].set_title(str(p))
    axs[1][count].hist(absav,bins=40, range=[0,1], stacked=True, fill=False)
    axs[1][count].set_title(str(p))
    if(count==2):
        title="Correlation " + "N=" + str(N) + " T="+str(T) + "(up)" + "        Averages (down)"
        m.suptitle(title)
    return axs

# Function to plot nodes neighbor list (set printnet=1)

def netprint(node_neigh_list):
    for i,w in enumerate(node_neigh_list):
        rede.write("{:d} - ".format(i))
        for j in w:
            rede.write("{:d} ".format(j))
        rede.write("\n")
               
# Here starts the main program       
global a,b,beta,c,dt,noise,N
a=-1
b=2.
c=-0.9
icto=math.sqrt((-b+math.sqrt(b**2-4*a*c))/(2*a))
beta=1.4
dt=0.01
noise=0.04
N=20
T=50
Nij=N*(N-1)/2
printhisto=1 #set to 1 to print the histogram
printnet=1 #set to 1 to print the nodes neighbor list
fig, axs = m.subplots(nrows=2,ncols=6)
saida=open("z6.out","w")
saida.write("#p, cor. media, desvio cor. media, ictogenicity \n")
rede=open("z6-net.out","w")
count=0
for p in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    print 'p=',p
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
        g = nx.erdos_renyi_graph(N,p)
    
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
    avicto=0.0
    cov=list(0 for i in range(Nij))
    avphi=list(0 for i in range(N))
    while(t<T):
        t+=dt
        if(t%10<dt):
            print "t=",int(t)
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
    if(printhisto==1):
        histograma(abscov,absav,count)
        count+=1

    avabscov=sum(abscov)/Nij
    desvio_correlacao=sum(map(lambda i:(i-avabscov)**2, abscov))/Nij
    saida.write('{:f} {:f} {:f} {:f}\n'.format(p,avabscov,desvio_correlacao,avicto)) 
#Graph
#    m.subplot(2,1,1)
#m.plot(x,b)
#    m.subplot(2,1,2)
#m.plot(x,abscov)
m.show()    


