import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from asdf import AsdfFile
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
import multiprocessing
import time
#from playsound import playsound
import math
def read_1(n):    
    files_path=[]
    flname = "{:0>{}}".format(n, 3)
    files_path.append(path+"halo_info_"+str(flname)+'.asdf')
    #cat=(CompaSOHaloCatalog(files_path, subsamples=dict(A=True)))
    #cat=CompaSOHaloCatalog(files_path, subsamples=dict(A=True, rv=True, pid=True), unpack_bits=True,cleaned=True, filter_func=lambda h: h['N'] >= 5000)
    cat=CompaSOHaloCatalog(files_path, subsamples=dict(A=True, rv=True, pid=True), unpack_bits=True,cleaned=True)
    hal={key: cat.halos[key] for key in ['id', 'npstartA','npoutA','x_com']}
    cat.subsamples = {key: cat.subsamples[key] for key in ['pid', 'lagr_pos', 'lagr_idx', 'tagged', 'density', 'pos', 'vel']}
    subs=cat.subsamples
    header=cat.header
    #cat.halos = {key: cat.halos[key] for key in ['id', 'npstartA','npoutA','x_com', 'x_L2com', 'v_com','N_merge', 'is_merged_to','haloindex','SO_central_density']}
    #return(np.column_stack((np.array(cat.subsamples['pos'][:,0]),np.array(cat.subsamples['pos'][:,1]),np.array(cat.subsamples['pos'][:,2]))))
    cat=0
    return([hal,subs,header])

def plot3d(x,y,z,N,a):
    s=time.time()
    z_ind=np.where((z >= -N) & (z <= N))
    y_ind=np.where((y >= -N) & (y <= N))
    x_ind=np.where((x >= -N) & (x <= N))
    idx_1=np.intersect1d(x_ind,y_ind)
    idx=np.intersect1d(idx_1,z_ind)
    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection ="3d")
    ax.scatter3D(x[idx],y[idx],z[idx],alpha=a)
    plt.title("Subsample Particles")
    e=time.time()
    print("time=%5.6f"%((e-s)*1000))
    
    
    
def plot3d_2(x,y,z,N,a):
    s=time.time()
    z_ind=np.where((z >= -N) & (z <= N))
    y_ind=np.where((y >= -N) & (y <= N))
    x_ind=np.where((x >= -N) & (x <= N))
    idx_1=np.intersect1d(x_ind,y_ind)
    plt.axes(projection ="3d")
    plt.scatter3D(x[idx],y[idx],z[idx],alpha=a)
    plt.title("Subsample Particles")
    e=time.time()
    print("time=%5.6f"%((e-s)*1000))
    

def slab_finder(x,box_size,num_of_file):
    def num_file(path):
        if (num_of_file==0):
            import os
            path2 =path
            count = 0
            for path in os.listdir(path2):
                if os.path.isfile(os.path.join(path2, path)):
                    count += 1
            return(count)
        else:
            return(num_of_file+1)
    N_f=int(num_file(path))-1
    N_per_f=box_size/N_f
    slab=int((box_size/2+x)/N_per_f)
    return(slab)

def plot_halo(halo_id):
    sub_pos=cat.subsamples['pos'][cat.halos['npstartA'][halo_id]:cat.halos['npstartA'][halo_id]+cat.halos['npoutA'][halo_id]]
    halo_pos=cat.halos["x_com"][halo_id]
    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection ="3d")
    ax.scatter3D(halo_pos[0],halo_pos[1],halo_pos[2],alpha=1,color='black')
    ax.scatter3D(sub_pos[:,0],sub_pos[:,1],sub_pos[:,2],alpha=.01,color='purple')
    
def read3(path,start,end):    
    files_path=[]
    for i in range(start,end+1):
        flname = "{:0>{}}".format(i, 3)
        files_path.append(path+"halo_info_"+str(flname)+'.asdf')
    cat=CompaSOHaloCatalog(files_path, subsamples=dict(A=True, rv=True, pid=True), unpack_bits=True,cleaned=True, filter_func=lambda h: h['N'] >= 5000)
    #cat=CompaSOHaloCatalog(files_path, subsamples=dict(A=True, rv=True, pid=True), unpack_bits=True,cleaned=True)
#    cat=CompaSOHaloCatalog(files_path)
#    x=cat.subsamples['pos'][:,0]
#    y=cat.subsamples['pos'][:,1]
#    z=cat.subsamples['pos'][:,2]
    #df = pd.DataFrame({"x" : x, "y" : y,"z":z})
    name="data"+str(start)+"_to_"+str(end)+".csv"
    return(cat)


# Estimate average density at a distance r from the CoM of halo
def density(r,r_value,s,cat_dens):
    denAt_r=np.array([0])
    dr=r[1]-r[0]
    for i in r:
        dens_ind = np.where((r_value >= i) & (r_value <= i+dr))
        dens_ind = list(s+dens_ind[0])
        #print(dens_ind)
        dens=cat_dens[dens_ind]
        #print(dens)
        if((len(dens)==0) & (i < max(r)/4)):
            denAt_r=np.append(denAt_r,denAt_r[len(denAt_r)-1])
            #denAt_r=np.append(denAt_r,0)
            continue
        l=0
        
        # Terminate the loop when a zero density region is found
        if len(dens)==0:
            l=np.where(r==i)
            l=l[0][0]
            break
        avg_dens=sum(dens)/len(dens)
        '''if avg_dens==0:
            l=np.where(r==i)
            l=l[0][0]
            break'''
        denAt_r=np.append(denAt_r,avg_dens)
    return(denAt_r,l)


def nfw(x,rho,r_s):
    return (rho/((x/r_s)*(1+x/r_s)**2))

#def nfw(x,rho,r_s):
#    return (rho/(((x/r_s)**1.5)(1+(x/r_s)**1.5)))
    
# Define the cost function to minimize
def cost_function(params):
    rho, r_s= params
    predicted_y = nfw(x,rho, r_s)
    mse = np.mean((predicted_y - y)**2)
    return mse
def kernel(distances):
    bandwidth = 0.1
    weights = np.exp(-0.5 * (distances / bandwidth)**2)
    return weights / np.sum(weights)

params2=np.array([0,0])      

'''def fit(r, denAt_r,l):
    #m=np.where(denAt_r==max(denAt_r))
    #m=m[0][0]
    #x=(r[m:l])
    #y=(denAt_r[m:l])
    # Calculate the distances between each data point and all others
    distances = np.abs(x[:, np.newaxis] - x[np.newaxis, :])

    # Calculate the weights using the kernel function
    weights = kernel(distances)

    # Estimate the parameters that minimize the cost function
    initial_params = [.1, .1]
    result = minimize(cost_function, initial_params)
    rho, r_s = result.x
    print(rho,r_s)
    #parms2=np.vstack((params2,np.array([rho,r_s])))
    # Evaluate the function with the estimated parameters
    y_fit = nfw(x, rho, r_s)
    return([x,y,y_fit,[rho,r_s]])'''


from scipy.optimize import curve_fit

def fit():

    popt, pcov = curve_fit(nfw, x, y, bounds=(0, [11000, 5.]))
    return(popt)
    
def merge(Af):
    class catalogue:
        def __init__(self,halos, subsamples, header):
            self.halos=halos
            self.subsamples=subsamples
            self.header=header
        
    halo_dict=Af[0][0]
    #Af[0][0]=0

    #Af[0][1]=0
    header=Af[0][2]
    subsample_dict=(Af[0][1])
    for k in range(1,len(Af)):
        #halo_dict= merg`eDictionary(halo_dict, Af[k][1])
        #subsample_dict=mergeDictionary(subsample_dict, Af[k][1])
        for key in halo_dict.keys():
            # Combine the dictionaries for each key
            halo_dict[key] = np.concatenate((np.array(halo_dict[key]),np.array(Af[k][0][key])), axis=0)  
        for key2 in subsample_dict.keys():
            #print(key2)
            subsample_dict[key2] = np.concatenate((np.array(subsample_dict[key2]),np.array(Af[k][1][key2])), axis=0)
    cat=catalogue(halo_dict,subsample_dict,header)
    
    return cat #[halo_dict,subsample_dict,header]

def read(start,end,read_1):
    start1 = time.time()
    pool=multiprocessing.Pool()
    inputs = np.arange(end-start)+start

    Af = pool.map(read_1, inputs)
    end1 = time.time()
    print("The time of execution:",(end1-start1) * 10**3, "ms")
    return(merge(Af))



def plot_halo(idx):
    
    s=cat.halos["npstartA"][idx]
    M=cat.halos["npoutA"][idx]
    x_h=cat.halos["x_com"][idx]
    
    x_b=cat.subsamples['pos'][s:s+M,0]
    y_b=cat.subsamples['pos'][s:s+M,1]
    z_b=cat.subsamples['pos'][s:s+M,2]
    
    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection ="3d")
    ax.scatter3D(x_b,y_b,z_b,alpha=.1)
    ax.scatter3D(x_h[0],x_h[1],x_h[2])
    plt.title("Subsample Particles")
    plt.show()
def vel(r,rho,r_s):
    G= 4.3009172706*10**(-9) #Mpc x Msun/(km/s)
    l= 1#10**(-6)
    #return(np.sqrt(4*np.pi*G*rho*r_s**3/r*(np.log(1+r/r_s)-r/r_s/(r/r_s+1))))
    return(np.sqrt(l*4*np.pi*G*rho*r_s**3/r*(np.log(1+r/r_s)-r/r_s/(r/r_s+1))))#*cat.header["H0"]/100)


def convert(x0,x): 
    return (np.array([np.sqrt(sum((x0-x)*(x0-x)))]))

def multi2(t,x0,k):
    A=np.array([convert(x0,B[t])])
    for i in range(t+1,t+k):
        A=np.append(A,np.array([convert(x0,B[i])]),axis=0)
    return(np.array(A))

def x2r(x0):
    import cupy as cp
    import math
    st=time.time()
    
    #print(np.array(B))
    nthread=16
    
    k=math.ceil(len(B)/nthread)-1
    #print(k)
    #rang=np.arange(nthread)*k
    rang= [(j*k, x0, k) for j in range(nthread)] # 0,3,....45
    #print(rang)
    with multiprocessing.Pool() as pool:
        #result = pool.map(multi2, rang)
        result = pool.starmap(multi2, rang)
    #print(result)
    #print(result)
    result=np.vstack(result)
    #print(f"took time {time.time()-st}")
    r_rem=np.array(multi2(k*nthread,x0,len(B)-k*nthread))
    r_value=np.concatenate((result,r_rem),axis=0)
    #r_value=np.array(r_value.get())
    r_value=np.vstack(r_value)
    #return (result)
    return(r_value)


from scipy.optimize import fsolve
import numpy as np

# Define the equation

def f_rs(rs, r1, r2, m1, m2):
    #r1, r2, m1, m2 = 2, 3, 10, 15
    return (math.log((rs + r2) / rs) - r2 / (rs + r2)) * m1 / m2 - (math.log((rs + r1) / rs) - r1 / (rs + r1))

def f_rho_s(rs, r1, m1):
  return m1 / (
      (
          math.log((rs + r1) / rs) - r1 / (rs + r1)
      ) * 4 * math.pi * rs**3
  )
def solve():
    solution = fsolve(equation, .1)
def plot_halos(mass_range,samp_frac):    
    import random
    import math
    #Index of the halo

    #samp_frac= 5  # percentage
    i=0
    col=3
    fig, ax3 = plt.subplots(nrows=(math.ceil(len(mass_range)/col)), ncols=col,figsize=(18, int(1.5*len(halos))),subplot_kw={'projection': '3d'})
    #fig = plt.figure(figsize=plt.figaspect(0.5))
    #ax3 = fig.add_subplot(projection='3d')

    for idx in mass_range:
        s=cat.halos["npstartA"][idx]
        M=cat.halos["npoutA"][idx]
        x_h=cat.halos["x_com"][idx]

    #    x_b=np.array(cat.subsamples['pos'][s:s+M,0])
    #    y_b=cat.subsamples['pos'][s:s+M,1]
    #    z_b=cat.subsamples['pos'][s:s+M,2]
        rand_pos=np.array(random.sample(list(cat.subsamples['pos'][s:s+M]), int(len(cat.subsamples["pos"][s:s+M])*samp_frac/100)))
        x_b=rand_pos[:,0]
        y_b=rand_pos[:,1]
        z_b=rand_pos[:,2]
        plt.title("Subsample Particles")
        #ax3.scatter3D(x_b,y_b,z_b,alpha=.1)
        ax3[int(i/col),int(i%col)].scatter(x_b,y_b,z_b,alpha=.1)
    #    ax3[int(i/col),int(i%col)].scatter3D(x_h[0],x_h[1],x_h[2])
    #    ax3[int(i/col),int(i%col)].set_xlabel("x")
    #    ax3[int(i/col),int(i%col)].set_ylabel("y")
    #    ax3[int(i/col),int(i%col)].set_zlabel("z")
        #ax2[int(i/col),int(i%col)].legend(loc='upper right')
    #    fig.suptitle('Velocity distribution')
        i=i+1
    plt.show()
