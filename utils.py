
import numpy as np
import argparse

from src import mashf90


class CustomArgumentParser(argparse.ArgumentParser):
    """ Allow commented lines in input file """
    def convert_arg_line_to_args(self, line):
        if line.strip().startswith('#'):
            return []
        return line.split()

""" Unit conversions to atomic units """
cmm1 = 4.556335e-6
fs = 41.34136e0
kB = 3.166829e-6
eV = 1./27.2113961

def read_args():
    parser = CustomArgumentParser(fromfile_prefix_chars=["@","+"])
    conv = parser.add_argument_group("convergence")
    debye = parser.add_argument_group("debye")
    tully = parser.add_argument_group("tully")
    parser.add_argument("-model", type=str, default="spinboson", help="Model system", 
                        choices=["spinboson","biexciton","fmo3","fmo7","pyrazine","tully1","tully2","lh2"])
    parser.add_argument("-basis", type=str, default="site", help="Diabatic basis for certain model systems", 
                        choices=["exc","site","adia"])
    parser.add_argument("-units", type=str, default="au",help="""Choose unit system. 
                        cmm1: input in cmm1, output in fs. 
                        fs: input in fs, output in fs. """,
                        choices=["au","cmm1","fs"])
    parser.add_argument("-obstyp", type=str, default="pop", help="Observable types", 
                        choices=["pop","all","mannrich","nuc"])
    parser.add_argument("-init", type=int,help="Initial state (in Python indexing)")
    parser.add_argument("-initbasis",type=str,default='dia',choices=["dia","adia","exc"],help='Basis for initial state')
    parser.add_argument("-nucsamp",type=str,default='cl', help="Nuclear sampling (in case 'init' is not None). GS=ground state, WP=wavepacket",
                        choices=['classical','cl','wigner','wig','GS','WP'])
    parser.add_argument("-elsamp",type=str,default='theta',choices=["focused","theta"],help='Choice of initial distribution')
    parser.add_argument("-debug",action="store_true",help='Toggle debug code (plot energy conservation etc)')
    parser.add_argument("-mode",type=str,default='mash', help="Choose which version of SH to run (fssh not fully supported)",
                        choices=['mash','fssh'])
    parser.add_argument("-beta",type=float,default=1.,help="Reciprocal temperature [in a.u.]")
    parser.add_argument("-T",type=float,default=0,help="Temperature [in kelvin]")
    conv.add_argument("-dt",default=41,type=float,help="Time step")
    conv.add_argument("-nt","-TS",type=int,help="Number of time steps")
    conv.add_argument("-ntraj","-traj",type=int,help="Number of trajectories")
    conv.add_argument("-npar","-n",default=1,type=int,help="Number of parallel processors")
    debye.add_argument("-Delta",type=float,default=1.,help="Diabatic coupling")
    debye.add_argument("-epsilon","-eps",type=float,default=0.,help="(Half) energy bias")
    debye.add_argument("-lamda",type=float,default=0,help="System-bath reorganization energy.")
    debye.add_argument("-omegac",type=float,default=0.,help="Cutoff frequency.")
    debye.add_argument("-nf",type=int,default=1,help="Number of bath modes.")
    tully.add_argument("-WPenergy",type=float,default=0.,help="Wavepacket energy")
    tully.add_argument("-gamma",type=float,default=0.,help="Wavepacket width parameter")
    tully.add_argument("-pinit",type=float,default=0.,help="Initial momentum")
    args = parser.parse_args()

    """ ======= Convert input arguments to atomic units ======="""
    if args.units=="cmm1":
        args.epsilon = args.epsilon * cmm1
        args.Delta = args.Delta * cmm1
        args.lamda = args.lamda * cmm1
        args.omegac = args.omegac * cmm1
        args.dt = args.dt * fs
    if args.units=='fs':
        args.dt = args.dt * fs
    if args.T:
        args.beta = 1./(kB*args.T)
    
    return args


"""======== Initialize system ========="""
def setup_model(args):
    model = args.model
    epsilon = args.epsilon
    Delta = args.Delta
    nf = args.nf
    omega = np.zeros(nf,dtype=np.float64)
    if model=='spinboson':
        ns = 2
        Vconst = np.array([[epsilon,Delta],[Delta,-epsilon]])
        Vlin = np.zeros((nf,ns,ns),dtype=np.float64)
        fac = np.sqrt(0.5*args.lamda/nf)
        c = np.zeros(nf)
        for i in range(nf):
            omega[i] = args.omegac*np.tan(0.5*np.pi*(i+0.5)/nf)
            c[i] = fac*omega[i]
            Vlin[i,0,0] = c[i]
            Vlin[i,1,1] = -c[i]
    elif model=='fmo3':
        ns = 3
        Vlin = np.zeros((nf,ns,ns),dtype=np.float64)
        omega = np.zeros(nf,dtype=np.float64)
        Vconst = np.diag([12410.,12530,12210])
        Vconst[0,1] = -87.7
        Vconst[0,2] =   5.5
        Vconst[1,2] =  30.8
        Vconst = Vconst + Vconst.T - np.diag(np.diag(Vconst))
        Vconst -= np.eye(ns)*np.min(np.diag(Vconst))

        Vconst = Vconst*cmm1
        nf_site = nf//ns
        fac = np.sqrt(2*args.lamda/nf_site)
        c = np.zeros(nf_site)
        for i in range(nf_site):
            omega[i] = args.omegac*np.tan(0.5*np.pi*(i+0.5)/nf_site)
            c[i] = fac*omega[i]
        for i in range(ns):
            omega[i*nf_site:(i+1)*nf_site] = omega[:nf_site]
            Vlin[i*nf_site:(i+1)*nf_site,i,i] = -c
    elif model=='fmo7':
        ns = 7
        Vlin = np.zeros((nf,ns,ns),dtype=np.float64)
        omega = np.zeros(nf,dtype=np.float64)
        Vconst = np.diag([12410.,12530,12210,12320,12480,12630,12440])
        Vconst[0,1] = -87.7
        Vconst[0,2] =   5.5
        Vconst[0,3] = - 5.9
        Vconst[0,4] =   6.7
        Vconst[0,5] = -13.7
        Vconst[0,6] = - 9.9
        Vconst[1,2] =  30.8
        Vconst[1,3] =   8.2
        Vconst[1,4] =   0.7
        Vconst[1,5] =  11.8
        Vconst[1,6] =   4.3
        Vconst[2,3] = -53.5
        Vconst[2,4] = - 2.2
        Vconst[2,5] =  -9.6
        Vconst[2,6] =   6.0
        Vconst[3,4] = -70.7
        Vconst[3,5] = -17.0
        Vconst[3,6] = -63.3
        Vconst[4,5] =  81.1
        Vconst[4,6] = - 1.3
        Vconst[5,6] =  39.7
        Vconst = Vconst + Vconst.T - np.diag(np.diag(Vconst))
        Vconst -= np.eye(ns)*np.min(np.diag(Vconst))
        Vconst = Vconst*cmm1
        nf_site = nf//ns
        fac = np.sqrt(2*args.lamda/nf_site)
        c = np.zeros(nf_site)
        for i in range(nf_site):
            omega[i] = args.omegac*np.tan(0.5*np.pi*(i+0.5)/nf_site)
            c[i] = fac*omega[i]
        for i in range(ns):
            omega[i*nf_site:(i+1)*nf_site] = omega[:nf_site]
            Vlin[i*nf_site:(i+1)*nf_site,i,i] = c

    """ ======= Initialize mass ======= """
    mass = np.ones(nf)

    """ ====== Initialize fortran potential ===== """
    if 'tully' in model:
        nf = 1
        ns = 2
        mass = np.ones(1)*2000. 
        mashf90.init_tully(model[-1],mass)
    elif model in ['fmo3','fmo7']:
        mashf90.init_frexc(mass,omega,Vconst,c,nf,nf_site,ns)  
    else:
        mashf90.init_linvib(mass,omega,Vconst,Vlin,nf,ns)

    return mass,omega,nf,ns

""" Debugging section to check energy conservation, plot adiabatic populations etc. """
def debug(args,mass,omega,nf,ns):
    beta = args.beta
    dt = args.dt
    nt = args.nt
    q0,p0,qe0,pe0 = sample(args,mass,omega,nf,ns)

    import matplotlib.pyplot as plt
    """ Energy conservation """
    q,p,qe,pe,Et,ierr=mashf90.runtrj(q0[:,0], p0[:,0], qe0[:,0], pe0[:,0], dt, nt, nf, ns)
    t = dt*np.arange(nt+1)
    Eref = nf/beta
    plt.plot(t,Et/Eref,'o',alpha=0.5)

    """ Adiabatic populations """
    plt.figure()
    qa = np.zeros_like(qe)
    pa = np.zeros_like(pe)
    popa = np.zeros_like(qe)
    for it in range(nt+1):
        qa[it],pa[it] = mashf90.dia2ad(q[it],qe[it],pe[it])
        popa[it] = qa[it]**2 + pa[it]**2
    for n in range(ns):
        label=str(n) 
        plt.plot(t,popa[:,n],'-',color='C%i'%n,alpha=0.5,label=label)
    plt.legend()

    plt.show()

def sample(args,mass,omega,nf,ns):
    """ Sample npar sets of initial phase-space variables.
        qe and pe are real and imaginary parts of c in the diabatic representation """
    npar = args.npar
    beta = args.beta
    q = np.empty((nf,npar),order="F")
    p = np.empty((nf,npar),order="F")
    qe = np.empty((ns,npar),order="F")
    pe = np.empty((ns,npar),order="F")
    for j in range(npar):
        """ Nuclear sampling """
        if args.nucsamp in ['cl','classical']:
            """ Classical """
            q[:,j] = np.random.normal(0,1./np.sqrt(beta*mass*omega**2),nf)
            p[:,j] = np.random.normal(0,np.sqrt(mass/beta),nf)
        elif args.nucsamp in ['wig','wigner']:
            """ Wigner """
            qsig = np.sqrt(1./(2*mass*omega*np.tanh(beta*omega/2)))
            psig = np.sqrt(mass*omega/(2*np.tanh(beta*omega/2)))
            q[:,j] = np.random.normal(0,qsig,nf)
            p[:,j] = np.random.normal(0,psig,nf)
        elif args.nucsamp=='GS':
            """ Ground state """
            qsig = np.sqrt(1./(2*mass*omega))
            psig = np.sqrt(mass*omega/2.)
            q[:,j] = np.random.normal(0,qsig,nf)
            p[:,j] = np.random.normal(0,psig,nf)
        elif args.nucsamp=='WP':
            """ Wave packet """
            if args.gamma>0:
                qsig = 1./np.sqrt(2.*args.gamma)
                psig = np.sqrt(args.gamma/2.)
            else:
                qsig = psig = 0
            qmean = -15.
            if args.WPenergy:
                pmean = np.sqrt(2.*mass*args.WPenergy)
            else:
                pmean = args.pinit
            q[:,j] = np.random.normal(qmean,qsig,nf)
            p[:,j] = np.random.normal(pmean,psig,nf)

        """ Electronic sampling """
        assert(args.init<ns and args.init>=0)
        if args.elsamp=='theta':
            """ Theta sampling: sample from region where Theta_init(c) = 1 """
            while 1:
                normals = np.random.randn(2*ns)
                normals = normals/np.linalg.norm(normals)
                qej=normals[:ns]
                pej=normals[ns:]
                pop = qej**2+pej**2
                l=0
                for k in range(ns):
                    if pop[k]>pop[l]: l=k # Find index of state with largest population
                if l==args.init:
                    qe[:,j] = qej; pe[:,j] = pej
                    break
        elif args.elsamp=='focused':
            """ Focussed sampling: sample from points where Phi_n = \delta_{n,init} """
            Hn = np.sum(1./np.arange(1,ns+1))
            A = (ns-1.)/(Hn-1.)
            Pn = np.ones(ns)*(A-1.)/(A*ns)
            Pn[args.init] = 1./ns + (ns-1.)/(A*ns)
            phi = np.random.uniform(0,2*np.pi,ns)
            qe[:,j] = np.cos(phi)*np.sqrt(Pn)
            pe[:,j] = np.sin(phi)*np.sqrt(Pn)

        """ Perform basis transformation (if requested) """
        if args.initbasis=='exc':
            """ Convert from exc to site basis: transform given by q=0 """
            vad,U = mashf90.get_vad(0.*q[:,j],ns)
            qe[:,j] = np.dot(U,qe[:,j])
            pe[:,j] = np.dot(U,pe[:,j])
        elif args.initbasis=='adia':
            """ Convert from adia to dia """
            vad,U = mashf90.get_vad(q[:,j],ns)
            qe[:,j] = np.dot(U,qe[:,j])
            pe[:,j] = np.dot(U,pe[:,j])
    return q,p,qe,pe

def savedata(B,t,ns,args):
    tout=t/fs if args.units in ['cmm1','fs'] else t
    if args.obstyp=='pop':
        out = np.column_stack([tout,B.reshape((-1,ns))])
    elif args.obstyp=='all':
        out = np.column_stack([tout,B.reshape((-1,ns**2))])
    np.savetxt('%s.out'%args.obstyp,out)

