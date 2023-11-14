#!/usr/bin/env python3

import numpy as np
import os
import sys

import utils
from src import mashf90

"""=========== Read input parameters=========""" 
args = utils.read_args()

""" ==== Set random seed (comment line if this is not wanted) ==== """
np.random.seed(42)


""" ======= Save commonly used arguments in their own variables ======"""
model = args.model
beta = args.beta
nf = args.nf
dt = args.dt
nt = args.nt
ntraj = args.ntraj
npar = args.npar
obstyp = args.obstyp

""" Setup paralellization """
if npar>1:
    os.environ["OMP_NUM_THREADS"] = str(npar)

""" ==== Setup time """
t = np.arange(nt+1)*dt

"""======== Initialize potential ========="""
mass, omega, nf, ns = utils.setup_model(args)

"""===== Initialize MASH Fortran module===="""
mashf90.init_mash()

""" Debugging section to plot energy conservation, plot adiabatic populations etc. """
if args.debug:
    utils.debug(args,mass,omega,nf,ns)
    sys.exit()

"""Initialize observables"""
if obstyp=='pop':
    Bt = np.zeros((nt+1,ns))
elif obstyp=='all':
    Bt = np.zeros((nt+1,ns,ns),dtype=np.complex128)
elif obstyp=='nuc':
    qs = []; ps = []; ws = []
    Tpop = np.zeros((nt+1,ns))
    Rpop = np.zeros((nt+1,ns))


"""======== Loop over trajectories ========"""
ndiscarded = 0
for itraj in range(ntraj//npar):
    q,p,qe,pe = utils.sample(args,mass,omega,nf,ns)

    reps = {'site':'d','exc':'e','adia':'a','dia':'d'}
    rep = reps[args.basis]
    if obstyp=='pop':
        """ Measure population dynamics """
        bt, Et, ierr = mashf90.runpar_poponly(q, p, qe, pe, rep, dt, nt, nf, ns, npar)
    elif obstyp=='all':
        """ Measure dynamics of populations and coherences """
        bt, Et, ierr = mashf90.runpar_all(q, p, qe, pe, rep, dt, nt, nf, ns, npar)    
    elif obstyp=='nuc':
        """ Measure final nuclear distribution (Tully) """
        ierr = np.zeros(npar)
        for j in range(npar):
            qt,pt,qet,pet,Et,ierr[j]=mashf90.runtrj(q[:,j], p[:,j], qe[:,j], pe[:,j], dt, nt, nf, ns)
            qs.append(qt[-1,0])
            ps.append(pt[-1,0])

            popt = np.array([mashf90.mash_pops(qt[it],qet[it],pet[it],'d',2) for it in range(nt+1)])
            Tpop += (qt>0)*popt
            Rpop += (qt<0)*popt
            
    """ Check for failed trajectories """
    if sum(ierr)>0:
        ndiscarded += np.sum(ierr>0)
    
    """ Save observables """
    if args.obstyp in ['pop','all']:
        Bt += bt

    """ Store temporary results after each 10 % of the number of trajectories"""
    if ntraj > 10:
        if (itraj+1)%(ntraj//(10*npar)) == 0:
            ctraj = itraj*npar+1 - ndiscarded
            print(ctraj+ndiscarded)
            if args.obstyp in ['pop','all']:
                utils.savedata(Bt/ctraj,t,ns,args)
            np.savetxt('log.out',np.array([ctraj,ntraj]),fmt='%i')

""" Log number of successful trajectories as well as requested number of trajectories """
print('ndiscarded',ndiscarded)
ctraj = itraj*npar+1 - ndiscarded
np.savetxt('log.out',np.array([ctraj,ntraj]),fmt='%i')
ntraj = ctraj

""" Store final results """
if args.obstyp in ['pop','all']:
    Bt /= ntraj
    utils.savedata(Bt,t,ns,args)
if args.obstyp=='nuc':
    bins = 200
    qhist,qbins = np.histogram(qs,bins,density=True)
    phist,pbins = np.histogram(ps,bins,density=True)
    qbins = 0.5*(qbins[:-1]+qbins[1:])
    pbins = 0.5*(pbins[:-1]+pbins[1:])
    np.savetxt('qhist.out',np.column_stack([qbins,qhist]))
    np.savetxt('phist.out',np.column_stack([pbins,phist]))

    Tpop /= ntraj
    Rpop /= ntraj
    np.savetxt('scatt.out',np.column_stack([t/utils.fs,Tpop,Rpop]))
