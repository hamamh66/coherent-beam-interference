"""Check first-order two-mode absorption formulas. Requires numpy.
Finite Fock space is sufficient: the tested states have at most one photon
per input mode. These are algebra checks, not detector simulations.
"""
import numpy as np
rng=np.random.default_rng(171026)
a=np.array([[0,1],[0,0]],complex)
a1=np.kron(a,np.eye(2));a2=np.kron(np.eye(2),a)
n1=a1.conj().T@a1;n2=a2.conj().T@a2
vac=np.array([1,0,0,0],complex)
for _ in range(200):
 d1,d2=rng.normal(size=2)+1j*rng.normal(size=2)
 c=rng.normal(size=2)+1j*rng.normal(size=2);c/=np.linalg.norm(c)
 D=d1*a1+d2*a2  # channel operator b-hat = kappa1 a1 + kappa2 a2 in the manuscript
 psi=np.array([0,c[1],c[0],0],complex)
 assert np.allclose(D@psi,(d1*c[0]+d2*c[1])*vac)
 rate=lambda v:np.vdot(D@v,D@v).real
 assert np.allclose(rate(psi),abs(d1*c[0]+d2*c[1])**2)
 two=np.array([0,0,0,1],complex)
 assert np.allclose(rate(two),abs(d1)**2+abs(d2)**2)
 dark=np.array([0,-d1,d2,0],complex);dark/=np.linalg.norm(dark)
 assert np.allclose(D@dark,0)
 assert np.allclose(np.vdot(dark,(n1+n2)@dark),1)
 # Expanded expression for an arbitrary pure state in this exact sector.
 v=rng.normal(size=4)+1j*rng.normal(size=4);v/=np.linalg.norm(v)
 expanded=(abs(d1)**2*np.vdot(v,n1@v)+abs(d2)**2*np.vdot(v,n2@v)
           +2*np.real(d1.conjugate()*d2*np.vdot(v,a1.conj().T@a2@v)))
 assert np.allclose(rate(v),expanded)
print('PASS: 200 randomized checks of the channel expansion, single-excitation interference,')
print('two-mode |1,1> rate, and non-vacuum dark states. No physical experiment performed.')
