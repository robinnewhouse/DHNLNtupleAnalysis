import os
import ROOT
import sys
sys.path.append('../../python/')
import helpers
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors



f_br = helpers.ReadBRjson(os.path.dirname(os.path.abspath(__file__)) + '/BranchingRatios_DifferentMixings_Gronau_lifetime.json')

mass = ["2", "2.5", "3", "4", "4.5", "5", "7.5","10","12.5", "15", "17.5","20"]
model = ["e_only", "mu_only", "IH", "NH"]
channel = ["mmv", "mev", "emv", "eev"]

br_mu_only_mmv = []
br_mu_only_mev = []
br_mu_only_eev = []
br_mu_only_emv = []

br_e_only_eev = []
br_e_only_emv = []
br_e_only_mmv = []
br_e_only_mev = []

br_IH_eev = []
br_IH_emv = []
br_IH_mmv = []
br_IH_mev = []

br_NH_eev = []
br_NH_emv = []
br_NH_mmv = []
br_NH_mev = []

purple = matplotlib.colors.to_rgb('#86259B')
blue = matplotlib.colors.to_rgb('#01C0E0')
pink = matplotlib.colors.to_rgb('#CC0063')
orange = matplotlib.colors.to_rgb('#FE9603')



for m in mass: 

    br_e_only_eev.append( f_br.get_BR("eee", m,"e_only", use_str = True) )
    br_e_only_emv.append( f_br.get_BR("eeu", m,"e_only", use_str = True) )
    br_e_only_mmv.append( f_br.get_BR("uuu", m,"e_only", use_str = True) )
    br_e_only_mev.append( f_br.get_BR("uue", m,"e_only", use_str = True) )

    br_mu_only_mmv.append( f_br.get_BR("uuu", m,"mu_only", use_str = True) )
    br_mu_only_mev.append( f_br.get_BR("uue", m,"mu_only", use_str = True) )
    br_mu_only_eev.append( f_br.get_BR("eee", m,"mu_only", use_str = True) )
    br_mu_only_emv.append( f_br.get_BR("eeu", m,"mu_only", use_str = True) )

    br_IH_eev.append( f_br.get_BR("eee", m,"IH", use_str = True) )
    br_IH_emv.append( f_br.get_BR("eeu", m,"IH", use_str = True) )
    br_IH_mmv.append( f_br.get_BR("uuu", m,"IH", use_str = True) )
    br_IH_mev.append( f_br.get_BR("uue", m,"IH", use_str = True) )

    br_NH_eev.append( f_br.get_BR("eee", m,"NH", use_str = True) )
    br_NH_emv.append( f_br.get_BR("eeu", m,"NH", use_str = True) )
    br_NH_mmv.append( f_br.get_BR("uuu", m,"NH", use_str = True) )
    br_NH_mev.append( f_br.get_BR("uue", m,"NH", use_str = True) )

plt.figure(0)
plt.plot(mass,br_e_only_eev, linestyle="", marker= 'o',label ="BR($N \\rightarrow ee\\nu$)",markersize=6, color = blue)
plt.plot(mass,br_e_only_emv, linestyle="", marker= 'o',label ="BR($N \\rightarrow e\\mu\\nu$ )",markersize=6, color = orange)
plt.plot(mass,br_e_only_mev, linestyle="", marker= 'o',label ="BR($N \\rightarrow \\mu e\\nu$)",markersize=6, color = pink)
plt.plot(mass,br_e_only_mmv, linestyle="", marker= 'o',label ="BR($N \\rightarrow \\mu\\mu\\nu$ )",markersize=6, color = purple)
plt.xlabel("HNL mass [GeV] ",fontsize=12)
plt.ylabel("BR($N \\rightarrow ll\\nu$) - any flavour neutrino",fontsize=12)
plt.title("Electron coupling only \n ($x_e$ = 1 , $x_\\mu$ = 0 , $x_\\tau$ = 0) ",fontsize=12)
plt.legend(numpoints=1,loc='upper right')

plt.savefig('br_e_only.pdf')

plt.figure(1)
plt.plot(mass,br_mu_only_eev, linestyle="", marker= 'o',label ="BR($N \\rightarrow ee\\nu$)",markersize=6, color = blue)
plt.plot(mass,br_mu_only_emv, linestyle="", marker= 'o',label ="BR($N \\rightarrow e\\mu\\nu$ )",markersize=6, color = orange)
plt.plot(mass,br_mu_only_mev, linestyle="", marker= 'o',label ="BR($N \\rightarrow \\mu e\\nu$)",markersize=6, color = pink)
plt.plot(mass,br_mu_only_mmv, linestyle="", marker= 'o',label ="BR($N \\rightarrow \\mu\\mu\\nu$ )",markersize=6, color = purple)
plt.xlabel("HNL mass [GeV] ",fontsize=12)
plt.ylabel("BR($N \\rightarrow ll\\nu$) - any flavour neutrino",fontsize=12)
plt.title("Muon coupling only \n ($x_e$ = 0 , $x_\\mu$ = 1 , $x_\\tau$ = 0)",fontsize=12)
plt.legend(numpoints=1,loc='upper right')

plt.savefig('br_mu_only.pdf')

plt.figure(2)
plt.plot(mass,br_IH_eev, linestyle="", marker= 'o',label ="BR($N \\rightarrow ee\\nu$)",markersize=6, color = blue)
plt.plot(mass,br_IH_emv, linestyle="", marker= 'o',label ="BR($N \\rightarrow e\\mu\\nu$ )",markersize=6, color = orange)
plt.plot(mass,br_IH_mev, linestyle="", marker= '^',label ="BR($N \\rightarrow \\mu e\\nu$)",markersize=4, color = pink)
plt.plot(mass,br_IH_mmv, linestyle="", marker= 'o',label ="BR($N \\rightarrow \\mu\\mu\\nu$ )",markersize=6, color = purple)
plt.xlabel("HNL mass [GeV] ",fontsize=12)
plt.ylabel("BR($N \\rightarrow ll\\nu$) - any flavour neutrino",fontsize=12)
plt.title("Inverted neutrino hierarchy model \n ($x_e$ = 1/3 , $x_\\mu$ = 1/3 , $x_\\tau$ = 1/3)",fontsize=12)
plt.legend(numpoints=1,loc='upper right')
plt.savefig('br_IH.pdf')

plt.figure(3)
plt.plot(mass,br_NH_eev, linestyle="", marker= 'o',label ="BR($N \\rightarrow ee\\nu$)",markersize=6, color = blue)
plt.plot(mass,br_NH_emv, linestyle="", marker= 'o',label ="BR($N \\rightarrow e\\mu\\nu$ )",markersize=6, color = orange)
plt.plot(mass,br_NH_mev, linestyle="", marker= 'o',label ="BR($N \\rightarrow \\mu e\\nu$)",markersize=6, color = pink)
plt.plot(mass,br_NH_mmv, linestyle="", marker= 'o',label ="BR($N \\rightarrow \\mu\\mu\\nu$ )",markersize=6, color = purple)
plt.xlabel("HNL mass [GeV] ",fontsize=12)
plt.ylabel("BR($N \\rightarrow ll\\nu$) - any flavour neutrino",fontsize=12)
plt.title("Normal neutrino hierarchy model \n ($x_e$ = 0.06 , $x_\\mu$ = 0.48 , $x_\\tau$ = 0.46)",fontsize=12)
plt.legend(numpoints=1,loc='upper right')
plt.savefig('br_NH.pdf')