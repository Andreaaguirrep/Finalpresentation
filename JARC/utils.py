from ocelot import *
from ocelot.cpbd.elements import *
import matplotlib.pyplot as plt

def trunc_s(p_array, smax):

    x = []
    y = []
    px = []
    py = []
    tau = []
    p = []
    q = []
    
    for i in range(len(p_array.rparticles[0])):
        if abs(p_array.rparticles[4][i]) < smax:
            x.append(p_array.rparticles[0][i])
            px.append(p_array.rparticles[1][i])
            y.append(p_array.rparticles[2][i])
            py.append(p_array.rparticles[3][i])
            tau.append(p_array.rparticles[4][i])
            p.append(p_array.rparticles[5][i])
            q.append(p_array.q_array[i])


    p_array_truncated=ParticleArray(n=len(x))
    p_array_truncated.E = p_array.E # GeV

    p_array_truncated.rparticles[0] = np.array(x)
    p_array_truncated.rparticles[1] = np.array(px)
    p_array_truncated.rparticles[2] = np.array(y)
    p_array_truncated.rparticles[3] = np.array(py)
    p_array_truncated.rparticles[4] = np.array(tau)
    p_array_truncated.rparticles[5] = np.array(p)
    p_array_truncated.q_array = np.array(q)

    return p_array_truncated

def show_parameter_evolution(tws_track):
    sigma_x = np.sqrt([tw.xx for tw in tws_track]) * 1.e6
    sigma_y = np.sqrt([tw.yy for tw in tws_track]) * 1.e6
    s = np.array([tw.s for tw in tws_track])
    E = np.array([tw.E for tw in tws_track])

    fig, axs = plt.subplots(2)
    fig.suptitle('Beam parameter evolution')


    p1,=axs[0].plot(s, sigma_x,'b-')
    axs[0].set_xlabel("b (m)")
    axs[0].set_ylabel(r"$\sigma_x$, ($\mu$ m)")
    
    ax2 = axs[0].twinx()
    p2,=ax2.plot(s, sigma_y,'r-')
    ax2.set_xlabel("s (m)")
    ax2.set_ylabel(r"$\sigma_y$, ($\mu$ m)")
    axs[0].legend([p1,p2],[r'$\sigma_x$',r'$\sigma_y$'])

    axs[1].plot(s,E)
    axs[1].set_xlabel("s (m)")
    axs[1].set_ylabel("E (GeV)")