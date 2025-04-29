from ocelot import *
from ocelot.cpbd.elements import *
from ocelot.gui.accelerator import *
from ocelot.cpbd.track import track_debug
import matplotlib.pyplot as plt
import numpy as np


#from rfglinacbte import *
#from Jarc import *
from linaccopy import *

#lat = MagneticLattice(lattice_list)
lat = MagneticLattice(LINAC)   



##### PA1RFGUN #####
# AX =-21.12   
# BX =24.36 
# AY =-21.16     
# BY =24.42      
# EMITX =1.2023504614117648e-06
# EMITY =1.2023504614117648e-06 
##### PA1RFGUN #####


##### PSECTB #####
# AX =.8330530584703799  
# BX =8.804855185360497    
# AY =.8104567519439174
# BY =8.82335409438248       
# EMITX =10.2023504614117648e-06 
# EMITY =10.2023504614117648e-06 
##### PSECTB #####

AX =.8241693577069001 
BX =13.14870543090723   
AY =.3785236088247831
BY =7.530918524683517   
EMITX =10.2023504614117648e-06 
EMITY =10.2023504614117648e-06 

energy = 2.5 # GeV  #why?
tw0 = Twiss()
tw0.alpha_x = AX
tw0.beta_x = BX

tw0.alpha_y = AY
tw0.beta_y = BY
tw0.E = energy
tw0.emit_xn = EMITX 
tw0.emit_yn = EMITY
tw0.emit_x = EMITX   / energy * 0.511e-3
tw0.emit_y = EMITY / energy * 0.511e-3

tws = twiss(lat,tw0)

plot_opt_func(lat, tws,legend=False, grid=False, top_plot=['Dx','Dy'])


plt.show()
