#!/usr/bin/env python3

import numpy as np
import os
import sys

import model
from src import mashf90

"""=========== Read input parameters =========""" 
args = model.read_args()

""" Save commonly used arguments in their own variables """
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

""" Setup time """
t = np.arange(nt+1)*dt

"""======== Initialize potential ========="""
mass, omega, nf, ns = model.setup_model(args)

"""===== Initialize MASH Fortran module ===="""
mashf90.init_mash(beta)

""" Debugging section to plot energy conservation, plot adiabatic populations etc. """
if args.debug:
    model.debug(args,mass,omega,nf,ns)
    sys.exit()

""" ===== Initialize observables ==== """
if obstyp=='pop':
    Bt = np.zeros((nt+1,ns))
elif obstyp=='nuc':
    qs = []; ps = []; ws = []
    Tpop = np.zeros((nt+1,ns))
    Rpop = np.zeros((nt+1,ns))

"""======== Loop over trajectories ========"""
ndiscarded = 0
for itraj in range(ntraj//npar):
    q0,p0,qe0,pe0 = model.sample(args,mass,omega,nf,ns)
    q = np.array(q0.copy(),order='F')
    p = np.array(p0.copy(),order='F')
    qe = np.array(qe0.copy(),order='F')
    pe = np.array(pe0.copy(),order='F')
    reps = {'site':'d','exc':'e','adia':'a','dia':'d'}
    if obstyp=='pop':
        """ Measure population dynamics """
        rep = reps[args.basis]
        bt, Et, ierr = mashf90.runpar(q, p, qe, pe, rep, dt, nt, nf, ns, npar)
    elif obstyp=='nuc':
        """ Measure final nuclear distribution (Tully) """
        ierr = np.zeros(npar)
        for j in range(npar):
            qt,pt,qet,pet,at,Et,ierr[j]=mashf90.runtrj(q[:,j], p[:,j], qe[:,j], pe[:,j], dt, nt, nf, ns)
            qs.append(qt[-1,0])
            ps.append(pt[-1,0])
            popt = np.array([mashf90.mash_pops(qt[it],qet[it],pet[it],'d') for it in range(nt+1)])
            Tpop += (qt>0)*popt
            Rpop += (qt<0)*popt
    """ Check for failed trajectories """
    if sum(ierr)>0:
        ndiscarded += np.sum(ierr>0)
        print('ndiscarded',ndiscarded)
    
    """ Save observables """
    if args.obstyp in ['pop']:
        Bt += bt

    """ Store temporary results after each 10 % of the number of trajectories"""
    if ntraj > 10:
        if (itraj+1)%(ntraj//(10*npar)) == 0:
            ctraj = (itraj+1)*npar - ndiscarded
            print(ctraj+ndiscarded)
            if args.obstyp in ['pop']:
                model.savedata(Bt/ctraj,t,args,args.obstyp)
            np.savetxt('log.out',np.array([ctraj,ntraj]),fmt='%i')

""" Log number of successful trajectories as well as requested number of trajectories """
ctraj = (itraj+1)*npar - ndiscarded
np.savetxt('log.out',np.array([ctraj,ntraj]),fmt='%i')
ntraj = ctraj

""" Store final results """
if args.obstyp in ['pop']:
    Bt /= ntraj
    model.savedata(Bt,t,args,args.obstyp)
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
    np.savetxt('scatt.out',np.column_stack([t/model.fs,Tpop,Rpop]))
