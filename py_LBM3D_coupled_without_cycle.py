import numpy as np
#import matplotlib.pyplot as plt
#from DoRTE2D import FEMDoRTE2D
#from phasefuncISO import FEMphasefuncISO
#from phasefuncO import FEMphasefuncO
#from phasefuncAni import FEMphasefuncAni1
#from Fgettra3d import Fgettra3d
#from feqtra3d import feqtra3d


# # Placeholder for the custom functions
def FEMDoRTE2D(Ntheta, Nfai):
    # Implement this function based on your MATLAB code
    # Initialize arrays
    angletheta = np.zeros((Ntheta, 1))
    anglefai = np.zeros((Nfai, 1))
    weighttheta = np.zeros((Ntheta, 1))
    weightfai = np.zeros((Nfai, 1))

    deltatheta = np.pi / Ntheta
    deltafai = 2 * np.pi / Nfai

    m = np.arange(1, Ntheta + 1)
    angletheta = (m - 1/2) * deltatheta

    n = np.arange(1, Nfai + 1)
    anglefai = (n - 1/2) * deltafai

    anglethetahalfp = m * deltatheta  # +
    anglethetahalfm = (m - 1) * deltatheta  # -

    anglefaihalfp = n * deltafai  # +
    anglefaihalfm = (n - 1) * deltafai  # -

    # Calculate weights
    weighttheta = np.cos(anglethetahalfm) - np.cos(anglethetahalfp)
    weightfai = anglefaihalfp - anglefaihalfm

    return angletheta, anglefai, weighttheta, weightfai


def FEMphasefuncISO(mu):
    # Implement this function based on your MATLAB code
    sfx = 1 + (0*mu)
    return sfx

def FEMphasefuncO(mu):
    # Implement this function based on your MATLAB code
    sfx = 0*mu
    return sfx

def FEMphasefuncAni1(mu):
    # Implement this function based on your MATLAB code
    sfx = 1 + mu
    return sfx

def Fgettra3d(e, w, NX, NY, NZ, Q, mu, yita, kesi, c, lamda, Lx, Fd):
    # Implement this function based on your MATLAB code
    # Initialize the array d with zeros
    d = np.zeros((NX+1, NY+1, NZ+1, Q))
    
    # Iterate over the range of k from 0 to 14 (equivalent to MATLAB 1 to 15)
    for k in range(15):
        # Calculate the factor
        factor = (1 + lamda * 3 * e[k, 0] * mu * Lx / c +
                      lamda * 3 * e[k, 1] * yita * Lx / c +
                      lamda * 3 * e[k, 2] * kesi * Lx / c)
        
        # Update d array at the kth position
        d[..., k] = w[k] * Fd * factor
    
    return d


def feqtra3d(I0, e, w, NX, NY, NZ, Q, mu, yita, kesi, cs, c, Lx):
    # Implement this function based on your MATLAB code
    a = np.zeros((NX+1, NY+1, NZ+1, Q))
    d = np.zeros((3, 3, Q))

    b = np.zeros((3, 3))
    b[0, 0] = (mu * Lx) ** 2 - cs ** 2
    b[0, 1] = mu * Lx * yita * Lx
    b[0, 2] = mu * Lx * kesi * Lx
    b[1, 0] = yita * Lx * mu * Lx
    b[1, 1] = (yita * Lx) ** 2 - cs ** 2
    b[1, 2] = yita * Lx * kesi * Lx
    b[2, 0] = kesi * Lx * mu * Lx
    b[2, 1] = kesi * Lx * yita * Lx
    b[2, 2] = (kesi * Lx) ** 2 - cs ** 2

    for k in range(15):
        d[0, 0, k] = 3 * (e[k, 0] ** 2) - 1
        d[0, 1, k] = 3 * e[k, 0] * e[k, 1]
        d[0, 2, k] = 3 * e[k, 0] * e[k, 2]
        d[1, 0, k] = 3 * e[k, 0] * e[k, 1]
        d[1, 1, k] = 3 * (e[k, 1] ** 2) - 1
        d[1, 2, k] = 3 * e[k, 1] * e[k, 2]
        d[2, 0, k] = 3 * e[k, 2] * e[k, 0]
        d[2, 1, k] = 3 * e[k, 2] * e[k, 1]
        d[2, 2, k] = 3 * (e[k, 2] ** 2) - 1

    for k in range(15):
        a[:, :, :, k] = w[k] * I0 * (1 + 3 * e[k, 0] * mu * Lx / c + 3 * e[k, 1] * yita * Lx / c + 3 * e[k, 2] * kesi * Lx / c + \
                                     b[0, 0] * d[0, 0, k] / (2 * cs ** 2) + b[0, 1] * d[0, 1, k] / (2 * cs ** 2) + b[0, 2] * d[0, 2, k] / (2 * cs ** 2) + \
                                     b[1, 0] * d[1, 0, k] / (2 * cs ** 2) + b[1, 1] * d[1, 1, k] / (2 * cs ** 2) + b[1, 2] * d[1, 2, k] / (2 * cs ** 2) + \
                                     b[2, 0] * d[2, 0, k] / (2 * cs ** 2) + b[2, 1] * d[2, 1, k] / (2 * cs ** 2) + b[2, 2] * d[2, 2, k] / (2 * cs ** 2))

    feq = a
    return feq

# Angular Discretization
Ntheta = 8  # select Ntheta to be even number
Nfai = 16  # select Nfai to be even number
Nangle = Ntheta * Nfai
angletheta, anglefai, weighttheta, weightfai = FEMDoRTE2D(Ntheta, Nfai)

# mu = np.sin(angletheta) * np.cos(anglefai)
mu = np.outer(np.sin(angletheta), np.cos(anglefai)).flatten('F')
# eta = np.sin(angletheta) * np.sin(anglefai)
eta = np.outer(np.sin(angletheta), np.sin(anglefai)).flatten('F')
# xi = np.cos(angletheta) * (1 + 0 * anglefai)
xi = np.outer(np.cos(angletheta), (1+0*anglefai)).flatten('F')
mu = mu.reshape(Nangle, 1)
eta = eta.reshape(Nangle, 1)
xi = xi.reshape(Nangle, 1)
angle = np.hstack((mu, eta, xi))
weight1 = np.outer(weighttheta, weightfai).flatten('F')
angleweight = weight1.reshape(Nangle, 1)
N = [Ntheta, Nfai]  # Identity of angular discretization

# Phase function
ScatterScheme = 1
if ScatterScheme == 0:
    Phasefunc = FEMphasefuncO  # 0: no scattering
elif ScatterScheme == 1:
    Phasefunc = FEMphasefuncISO  # 1: isotropically
elif ScatterScheme == 2:
    Phasefunc = phasefuncAni1  # 1 + mu*mu: Anisotropically sfx = 1 + mu

NX, NY, NZ = 20, 20, 20

# Lattice coefficient D3Q15
Q = 15
e = np.array([
    [0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], 
    [0, 0, 1], [0, 0, -1], [1, 1, 1], [-1, -1, -1], [1, 1, -1], 
    [-1, -1, 1], [1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, 1, 1]
])
w = np.array([2/9] + [1/9] * 6 + [1/72] * 8)

Nsol = (NX + 1) * (NY + 1) * (NZ + 1)
dx = 1 / NX
Lx = dx * NX
Total = 5
dt = 0.5
maxT = Total / dt  # time step
c = dx / dt
cs = c / np.sqrt(3)
xelem = np.linspace(0, 1, NX + 1)

# Basic settings
nref = 1
Tref = 600
Stefanboltzmanzeta = 5.67e-8
Tao = 1  # optic length
Omega = 0.9  # albedo
beta = Tao / Lx
Ks = Omega * beta
Ka = (1 - Omega) * beta
Emall = 1
Es = 1
tau_r = 1
lamda = (tau_r - 0.5) / tau_r
N = 0.01
lamdac = (N * 4 * Stefanboltzmanzeta * Tref**3) / beta
density = 1e-2
cv = 1e-3
cl = 3e8
k1 = lamdac / (density * cv * Lx * cl)
tau_c = k1 / (dt * cs**2) + 1/2

# Initialization
T = 0.5 * Tref * np.ones((NX + 1, NY + 1, NZ + 1))
T[:, :, 0] = Tref
Ib = Stefanboltzmanzeta * nref**2 * T**4 / np.pi
I = np.zeros((NX + 1, NY + 1, NZ + 1, Nangle))
II = np.zeros((Nsol, Nangle))
feq = np.zeros((NX + 1, NY + 1, NZ + 1, Q, Nangle))
geq = np.zeros((NX + 1, NY + 1, NZ + 1, Q))
f = np.zeros((NX + 1, NY + 1, NZ + 1, Q, Nangle))
g = np.zeros((NX + 1, NY + 1, NZ + 1, Q))
F0 = np.zeros((NX + 1, NY + 1, NZ + 1, Q, Nangle))
F = np.zeros((NX + 1, NY + 1, NZ + 1, Q, Nangle))
Fd = np.zeros((NX + 1, NY + 1, NZ + 1, Nangle))
Scatter = np.zeros((NX + 1, NY + 1, NZ + 1, Nangle))
fcol = np.zeros((NX + 1, NY + 1, NZ + 1, Q, Nangle))
gcol = np.zeros((NX + 1, NY + 1, NZ + 1, Q))
G0 = np.zeros((NX + 1, NY + 1, NZ + 1))
T0 = np.zeros((NX + 1, NY + 1, NZ + 1))
Gn = II @ angleweight
G = Gn.reshape(NX + 1, NY + 1, NZ + 1)

# Initialize the radiative field
for cycle_r in range(1, 100):
    print(f'cycle_r = {cycle_r}')
    for anglei in range(Nangle):
        direction = angle[anglei, :]  
        mui, etai, xii = direction
        sidots = angle @ direction
        sfx = Phasefunc(sidots).reshape(sidots.size,1)

        Scatering = Ks / (4 / np.pi) * II @ (angleweight * sfx)
        #Scatter=reshape(Scatering,NX+1,NY+1,NZ+1);
        Scatter = Scatering.reshape(NX + 1, NY + 1, NZ + 1)
        Fd = -Lx * beta * I[:, :, :, anglei] + Lx * Scatter + Lx * Ka * Ib

        F[:, :, :, :, anglei] = Fgettra3d(e, w, NX, NY, NZ, Q, mui, etai, xii, c, lamda, Lx, Fd)

        for k in range(Q):
            fcol[:, :, :, k, anglei] = f[:, :, :, k, anglei] + (feq[:, :, :, k, anglei] - f[:, :, :, k, anglei]) / tau_r + dt * F[:, :, :, k, anglei] + dt * (F[:, :, :, k, anglei] - F0[:, :, :, k, anglei]) / 2

        for k in range(Q):
            f[:, :, :, k, anglei] = np.roll(fcol[:, :, :, k, anglei], e[k], axis=(0, 1, 2))

        I[:, :, :, anglei] = np.sum(f[:, :, :, :, anglei], axis=3)

        # Boundary treatment for radiative field
        I[0, :, :, anglei] = 2 * I[1, :, :, anglei] - I[2, :, :, anglei]
        I[NX, :, :, anglei] = 2 * I[NX - 1, :, :, anglei] - I[NX - 2, :, :, anglei]
        I[:, 0, :, anglei] = 2 * I[:, 1, :, anglei] - I[:, 2, :, anglei]
        I[:, NY, :, anglei] = 2 * I[:, NY - 1, :, anglei] - I[:, NY - 2, :, anglei]
        I[:, :, 0, anglei] = 2 * I[:, :, 1, anglei] - I[:, :, 2, anglei]
        I[:, :, NZ, anglei] = 2 * I[:, :, NZ - 1, anglei] - I[:, :, NZ - 2, anglei]

        if mui<0: #oyz
            I[NX-1, :, :, anglei]=Emall*Ib[NX-1,:,:]
        if mui>0: #oyz
            I[1, :, :, anglei]=Emall*Ib[1,:,:]

        if etai<0: #oxz
            I[:, NY-1, :, anglei]=Emall*Ib[:,NY-1,:]
        if etai>0: #oxz
            I[:, 1, :, anglei]=Emall*Ib[:,1,:]
        if xii<0: #oxy
            I[:, :, NZ-1, anglei]=Emall * Ib[:,:,NZ-1]
        if xii>0: #oxy
            I[:, :, 1, anglei]=Emall * Ib[:,:,1]

    F0=F

    for anglei in range(Nangle):
        RI=I[:,:,:,anglei].reshape(Nsol)
        II[:,anglei]=RI[:]
        direction=angle[anglei,:]; 
        mui, etai, xii = direction
        feq[:,:,:,:,anglei]=feqtra3d(I[:,:,:,anglei],e,w,NX,NY,NZ,Q,mui,etai,xii,cs,c,Lx)
        for k in range(Q):
            f[1,:,:,k,anglei]=f[2,:,:,k,anglei]+feq[1,:,:,k,anglei]-feq[2,:,:,k,anglei]
            f[NX-1,:,:,k,anglei]=f[NX,:,:,k,anglei]+feq[NX-1,:,:,k,anglei]-feq[NX,:,:,k,anglei]
            f[:,1,:,k,anglei]=f[:,2,:,k,anglei]+feq[:,1,:,k,anglei]-feq[:,2,:,k,anglei]
            f[:,NY-1,:,k,anglei]=f[:,NY,:,k,anglei]+feq[:,NY-1,:,k,anglei]-feq[:,NY,:,k,anglei]
            f[:,:,1,k,anglei]=f[:,:,2,k,anglei]+feq[:,:,1,k,anglei]-feq[:,:,2,k,anglei]
            f[:,:,NZ-1,k,anglei]=f[:,:,NZ,k,anglei]+feq[:,:,NZ-1,k,anglei]-feq[:,:,NZ,k,anglei]
    
    #convergen or not
    if (cycle_r% 10 ==0):
        # Gn = II * angleweight
        Gall = Gn.reshape(NX + 1, NY + 1, NZ + 1)
        # error = np.sum(np.abs(Gall - G0), axis=(0, 1, 2)) / np.sum(np.abs(Gall), axis=(0, 1, 2))
        # Calculate the numerator and denominator
        numerator = np.sum(np.abs(Gall - G0), axis=(0, 1, 2))
        denominator = np.sum(np.abs(Gall), axis=(0, 1, 2))
        print("numerator",numerator)
        print("denominator",denominator)

        # Handle the division, avoiding invalid values
        error = np.where(denominator == 0, np.inf, numerator / denominator)
        G0 = Gall
        if error<1e-6:
            break

        # II[:, anglei] = I[:, :, :, anglei].reshape(Nsol)

for i in range(NX+1):
    for j in range(NY+1):
        for ij in range(NZ+1):
            for k in range(Q):
                geq[i,j,ij,k]=w[k]*T[i,j,ij]
                g[i,j,ij,k]=geq[i,j,ij,k]


# begin the transient evolution process
for cycle in range(1, 100):
    print(f'cycle = {cycle}')
    k2=1/(density*cv*cl)
    #print(haji.shape)
    for k in range(Q):
        gcol[:, :, :, k] = g[:, :, :, k] + (geq[:, :, :, k] - g[:, :, :, k]) / tau_c - dt * w[k] * beta * (1-Omega) * k2 * (4*(T**4)*Stefanboltzmanzeta-G)
        #gcol(:,:,:,k)=g(:,:,:,k)+(geq(:,:,:,k)-g(:,:,:,k))./tau_c - dt.*w(k).*beta.*(1-Omega).*k2.*(4*T.^4*Stefanboltzmanzeta-G);
    for k in range(Q):
        g[:, :, :, k] = np.roll(gcol[:, :, :, k], e[k,:], axis=(0, 1, 2))

    #temperature boundary treatment
    T[1,:,:]=0.5*Tref
    T[NX-1,:,:]=0.5*Tref
    T[:,NY-1,:]=0.5*Tref
    T[:,1,:]=0.5*Tref
    T[:,:,NZ-1]=0.5*Tref
    T[:,:,1]=Tref

    #nonequilibrium extrapolation boundary scheme for temperature distribution function
    for i in range(NX+1):
        for j in range(NY+1):
            for ij in range(NZ+1):
                for k in range(Q):
                    geq[i,j,ij,k]=w[k]*T[i,j,ij]

    for k in range(Q):
        g[0, :, :, k] = g[1, :, :, k] + geq[0, :, :, k] - geq[1, :, :, k]
        g[NX, :, :, k] = g[NX - 1, :, :, k] + geq[NX, :, :, k] - geq[NX - 1, :, :, k]
        g[:, 0, :, k] = g[:, 1, :, k] + geq[:, 0, :, k] - geq[:, 1, :, k]
        g[:, NY, :, k] = g[:, NY - 1, :, k] + geq[:, NY, :, k] - geq[:, NY - 1, :, k]
        g[:, :, 0, k] = g[:, :, 1, k] + geq[:, :, 0, k] - geq[:, :, 1, k]
        g[:, :, NZ, k] = g[:, :, NZ - 1, k] + geq[:, :, NZ, k] - geq[:, :, NZ - 1, k]
    Ib = Stefanboltzmanzeta * nref**2 * T**4 / np.pi

    for anglei in range(Nangle):
        direction = angle[anglei, :]
        mui, etai, xii = direction
        sidots = angle @ direction
        sfx = Phasefunc(sidots).reshape(sidots.size,1)

        Scatering = Ks / (4 / np.pi) * II @ (angleweight * sfx)
        #Scatter=reshape(Scatering,NX+1,NY+1,NZ+1);
        Scatter = Scatering.reshape(NX + 1, NY + 1, NZ + 1)
        Fd = -Lx * beta * I[:, :, :, anglei] + Lx * Scatter + Lx * Ka * Ib

        F[:, :, :, :, anglei] = Fgettra3d(e, w, NX, NY, NZ, Q, mui, etai, xii, c, lamda, Lx, Fd)

        for k in range(Q):
            fcol[:, :, :, k, anglei] = f[:, :, :, k, anglei] + (feq[:, :, :, k, anglei] - f[:, :, :, k, anglei]) / tau_r + dt * F[:, :, :, k, anglei] + dt * (F[:, :, :, k, anglei] - F0[:, :, :, k, anglei]) / 2
            

        for k in range(Q):
            f[:, :, :, k, anglei] = np.roll(fcol[:, :, :, k, anglei], e[k], axis=(0, 1, 2))

        # boundary
        I[:, :, :, anglei] = np.sum(f[:, :, :, :, anglei], axis=3)

        # Boundary treatment for radiative field
        I[0, :, :, anglei] = 2 * I[1, :, :, anglei] - I[2, :, :, anglei]
        I[NX, :, :, anglei] = 2 * I[NX - 1, :, :, anglei] - I[NX - 2, :, :, anglei]
        I[:, 0, :, anglei] = 2 * I[:, 1, :, anglei] - I[:, 2, :, anglei]
        I[:, NY, :, anglei] = 2 * I[:, NY - 1, :, anglei] - I[:, NY - 2, :, anglei]
        I[:, :, 0, anglei] = 2 * I[:, :, 1, anglei] - I[:, :, 2, anglei]
        I[:, :, NZ, anglei] = 2 * I[:, :, NZ - 1, anglei] - I[:, :, NZ - 2, anglei]

        if mui<0: #oyz
            I[NX-1, :, :, anglei]=Emall*Ib[NX-1,:,:]
        if mui>0: #oyz
            I[1, :, :, anglei]=Emall*Ib[1,:,:]

        if etai<0: #oxz
            I[:, NY-1, :, anglei]=Emall*Ib[:,NY-1,:]
        if etai>0: #oxz
            I[:, 1, :, anglei]=Emall*Ib[:,1,:]
        
        if xii<0: #oxy
            I[:, :, NZ-1, anglei]=Emall*Ib[:,:,NZ-1]
        if xii>0: #oxy
            I[:, :, 1, anglei]=Emall*Ib[:,:,1]
    F0=F

    for anglei in range(Nangle):
        RI=I[:,:,:,anglei].reshape(Nsol)        
        II[:,anglei]=RI[:]
        direction=angle[anglei,:]; 
        mui, etai, xii = direction
        feq[:,:,:,:,anglei]=feqtra3d(I[:,:,:,anglei],e,w,NX,NY,NZ,Q,mui,etai,xii,cs,c,Lx)
        for k in range(Q):
            f[1,:,:,k,anglei]=f[2,:,:,k,anglei]+feq[1,:,:,k,anglei]-feq[2,:,:,k,anglei]
            f[NX-1,:,:,k,anglei]=f[NX,:,:,k,anglei]+feq[NX-1,:,:,k,anglei]-feq[NX,:,:,k,anglei]
            f[:,1,:,k,anglei]=f[:,2,:,k,anglei]+feq[:,1,:,k,anglei]-feq[:,2,:,k,anglei]
            f[:,NY-1,:,k,anglei]=f[:,NY,:,k,anglei]+feq[:,NY-1,:,k,anglei]-feq[:,NY,:,k,anglei]
            f[:,:,1,k,anglei]=f[:,:,2,k,anglei]+feq[:,:,1,k,anglei]-feq[:,:,2,k,anglei]
            f[:,:,NZ-1,k,anglei]=f[:,:,NZ,k,anglei]+feq[:,:,NZ-1,k,anglei]-feq[:,:,NZ,k,anglei]
        

    # Convergence judgment
    Gn = II @ angleweight
    G = Gn.reshape(NX + 1, NY + 1, NZ + 1)
    if (cycle_r% 10 ==0):
        #Gn = II * angleweight
        Gall = Gn.reshape(NX + 1, NY + 1, NZ + 1)
        error = np.sum(np.abs(Gall - G0), axis=(0, 1, 2)) / np.sum(np.abs(Gall), axis=(0, 1, 2))
        G0 = Gall
        if error<1e-6:
            break


import matplotlib.pyplot as plt
###############################################################
Rerr = np.linalg.norm(G - G0) / np.linalg.norm(G0)
print(f'Rerr = {Rerr}')
G0 = G
# if Rerr < 1e-6:
#     break

# Visualizing the results
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y, Z = np.meshgrid(np.linspace(0, 1, NX + 1), np.linspace(0, 1, NY + 1), np.linspace(0, 1, NZ + 1))

img = ax.scatter(X, Y, Z, c=G)
#img = ax.plot_surface(X,Y,Z, c=G)
fig.colorbar(img)
plt.show()
# # Save the plot to a file
# plt.savefig('3d_plot.png')
# plt.close()
##############################################################