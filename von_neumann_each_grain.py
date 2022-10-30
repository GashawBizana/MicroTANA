fig = plt.figure()
ax = fig.add_subplot()
dvdt=[]
gbn=[]

List=[i for i in range(79) if (i !=7 and i<12) or i >90]
# List=[5]
for kk in List:
    
    G_tc, G_dv, G_TL_l, G_TL_n, G_GB_n, G_mis, G_v, G_a, G_id, G_tcf,G_step,G_dis,G_tmc,G_af,G_tcf_p= [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    
    g=[kk,kk]
    tc=[-1000000000000,100000000000]
    dis=[-15000000,1500000]
    LS=0
    
    LE=len(Grain_VonNeuman_Data_for_timesteps)
    for i in range(0,LE):
        G_id0=[j[0] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]] 
        G_tc0=[j[1] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_dv0=[j[2] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_TL_l0=[j[3] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_TL_n0=[j[4] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_GB_n0=[j[5] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_mis0=[j[6] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_v0=[j[7] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_a0=[j[8] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_tcf0=[j[9] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_step0=[i for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_dis0=[j[12] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_tmc0=[j[13] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_af0=[j[10] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        G_tcf_p0=[j[15] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1]]
        
        G_id=G_id+G_id0
        G_tc=G_tc+G_tc0
        G_dv=G_dv+G_dv0
        G_TL_l=G_TL_l+G_TL_l0
        G_TL_n=G_TL_n+G_TL_n0
        G_GB_n=G_GB_n+G_GB_n0
        G_mis=G_mis+G_mis0
        G_v=G_v+G_v0
        G_a=G_a+G_a0
        G_tcf=G_tcf+G_tcf0
        G_step=G_step+G_step0
        G_dis=G_dis+G_dis0
        G_tmc=G_tmc+G_tmc0
        G_af=G_af+G_af0
        G_tcf_p=G_tcf_p+G_tcf_p0
        
    xb=-2*np.pi*np.array( G_tcf)*0.1/(np.array(G_af)*0.01)
    view=np.array(G_GB_n)
    dt=1/(step*0.01)
    xa1=np.array(G_tc)
    xa2=np.array(G_TL_l)
    xa=-(xa1-xa2/6)*0.1
    # xa=xa/(np.array(G_a) *0.01)
    # view=xa
    ya=np.array(G_dv)*0.001*dt
    ya=np.array(G_dis)*0.1*dt
    # xa=10*np.array(G_tcf_p)/np.array(G_af)
    
    ya=np.array(G_dis)*0.1
    # ya=np.array(G_v)*0.001
    xa=np.array(G_step)*0.01
    
    # xa=10*np.array(G_tc-xa2/6)/np.array(G_a)
    # print(np.mean(G_GB_n))
    # xa=-(np.array(G_tc)-xa2/6)
    # ya=np.array(G_dv)/np.array(G_af)
    # xa=-10*np.array(G_tc-xa2/6)/np.array(G_a)
    # xa=np.power(np.array(G_v)*0.001,-1/3)
    
    # xa=-(np.array(G_tcf))*0.1/np.array(G_af)
    # xa=np.array(G_TL_l)/6
    # plt.rcParams['text.usetex'] = True
    cmap=plt.cm.jet
    cmaplist=[cmap(i) for i in range(cmap.N)]
    cmaplist[0]=(0.5,0.5,0.5,1.0)
    cmap=mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist,cmap.N)
    
    # bounds=np.linspace(np.rint(min(mis)),np.rint(max(mis)),2*int(np.rint(max(mis))-np.rint(min(mis))))
    # bounds=np.linspace(min(view),max(view))
    bounds=np.linspace(0,38,37)
    
    # bounds=np.linspace(4,35)
    norm=mpl.colors.BoundaryNorm(bounds, cmap.N)
    # coef = np.polyfit(xa,ya,1)
    # poly1d_fn = np.poly1d(coef)
    # xnew=np.linspace(min(xa),max(xa),10)
    # fig = plt.figure()
    # ax = fig.add_subplot()
    
    # dvdt.append(coef[0])
    # gbn.append(G_GB_n[0])
    
    # sc=ax.scatter(xa,ya,s=50,marker='o',c=view,cmap=cmap,norm=norm)
    # ax.plot(xa,ya,alpha=1)
    ma = np.array(["o", "s", "d", "*",'+','>','<'])
    # if kk!=7:
    #     sc=ax.scatter(xa,np.cumsum(ya)-ya[0],marker=np.random.choice(ma),c=view,cmap=cmap,norm=norm)
    #     ax.plot(xa,np.cumsum(ya)-ya[0],alpha=0.8)
    
    # yy=np.cumsum(ya)-ya[0]
    yy=ya
    p = np.polyfit(xa,yy,2)
    p1 = np.polyfit(xa,xb,2)
    dp=np.polyder(p)
        
    sc=ax.scatter(np.polyval(p1, xa),np.polyval(p, xa),marker=np.random.choice(ma),c=view,cmap=cmap,norm=norm)
    ax.plot(np.polyval(p1, xa),np.polyval(p, xa),alpha=1,c='k')
    
    # sc=ax.scatter(xa,yy,marker=np.random.choice(ma),c=view,cmap=cmap,norm=norm)
    # ax.plot(xa,np.polyval(p, xa),alpha=1,c='k')
    # ax.plot(xa,view*100)
    
    # ax.set_xlabel('M/6-L ($nm$)')
    # ax.set_xlabel('k')
    # ax.set_ylabel('dV/dt ($nm^3$/step)')
    # ax.set_ylabel('V')
    # cb=plt.colorbar(sc)
    # cb.set_label('curvature')
    # cb.set_label('Number Of Faces')
    # plt.plot(xnew,poly1d_fn(xnew))
    # plt.axhline(y=0, color='r', linestyle='-',linewidth=2)
    # plt.axvline(x=0, color='k', linestyle='-',linewidth=2)
    # plt.axvline(x=-coef[1]/coef[0], color='k', linestyle='-',linewidth=2)
    # ax.scatter(xa,np.array(G_tc))
    
ax.set_ylabel('Relative position(nm)',fontweight='bold',size=12)
# ax.set_ylabel('Volume(nm^{3})',fontweight='bold',size=12)
ax.set_xlabel('t(ns)',fontweight='bold',size=12)
    
cb=plt.colorbar(sc)
# cb.set_label('Curvature (1/nm)',fontweight='bold',size=12)  
cb.set_label('Topological Class (n)',fontweight='bold',size=12)  
plt.show()


