import os
import ROOT
import sys
sys.path.append('../../python/')
import helpers
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors



f_br = helpers.ReadJsonFiles(os.path.dirname(os.path.abspath(__file__)) + '/BranchingRatios_DifferentMixings_Olegs_lifetime.json')

mass = ["2.0", "2.5", "3.0", "4.0", "4.5", "5.0", "7.5","10.0","12.5", "15.0", "17.5","20.0"]
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


mass_float= []
for m in mass: 
    mass_float.append(float(m))
    br_e_only_eev.append(f_br.get_br("eee", m, "e_only", use_str = True))
    br_e_only_emv.append(f_br.get_br("eeu", m, "e_only", use_str = True))
    br_e_only_mmv.append(f_br.get_br("uuu", m, "e_only", use_str = True))
    br_e_only_mev.append(f_br.get_br("uue", m, "e_only", use_str = True))

    br_mu_only_mmv.append(f_br.get_br("uuu", m, "mu_only", use_str = True))
    br_mu_only_mev.append(f_br.get_br("uue", m, "mu_only", use_str = True))
    br_mu_only_eev.append(f_br.get_br("eee", m, "mu_only", use_str = True))
    br_mu_only_emv.append(f_br.get_br("eeu", m, "mu_only", use_str = True))

    br_IH_eev.append(f_br.get_br("eee", m, "IH", use_str = True))
    br_IH_emv.append(f_br.get_br("eeu", m, "IH", use_str = True))
    br_IH_mmv.append(f_br.get_br("uuu", m, "IH", use_str = True))
    br_IH_mev.append(f_br.get_br("uue", m, "IH", use_str = True))

    br_NH_eev.append(f_br.get_br("eee", m, "NH", use_str = True))
    br_NH_emv.append(f_br.get_br("eeu", m, "NH", use_str = True))
    br_NH_mmv.append(f_br.get_br("uuu", m, "NH", use_str = True))
    br_NH_mev.append(f_br.get_br("uue", m, "NH", use_str = True))

plt.figure(0)
plt.plot(mass_float,br_e_only_eev, marker= 'o',label ="BR($N \\rightarrow ee\\nu$)",markersize=6, color = blue)
plt.plot(mass_float,br_e_only_emv, marker= 'v',label ="BR($N \\rightarrow e\\mu\\nu$ )",markersize=6, color = orange)
plt.plot(mass_float,br_e_only_mev, marker= 's',label ="BR($N \\rightarrow \\mu e\\nu$)",markersize=6, color = pink)
plt.plot(mass_float,br_e_only_mmv, marker= '^',label ="BR($N \\rightarrow \\mu\\mu\\nu$ )",markersize=6, color = purple)
plt.xlabel("$m_N$ [GeV] ",fontsize=14)
plt.ylabel("BR($N \\rightarrow ll\\nu$) - any flavour neutrino",fontsize=14)
plt.title("Electron-only mixing \n ($x_e$ = 1 , $x_\\mu$ = 0 , $x_\\tau$ = 0) ",fontsize=16)
plt.legend(numpoints=1,loc='upper right')

plt.savefig('br_e_only.pdf')

plt.figure(1)
plt.plot(mass_float,br_mu_only_eev, marker= 'o',label ="BR($N \\rightarrow ee\\nu$)",markersize=6, color = blue)
plt.plot(mass_float,br_mu_only_emv, marker= 'v',label ="BR($N \\rightarrow e\\mu\\nu$ )",markersize=6, color = orange)
plt.plot(mass_float,br_mu_only_mev, marker= 's',label ="BR($N \\rightarrow \\mu e\\nu$)",markersize=6, color = pink)
plt.plot(mass_float,br_mu_only_mmv, marker= '^',label ="BR($N \\rightarrow \\mu\\mu\\nu$ )",markersize=6, color = purple)
plt.xlabel("$m_N$ [GeV] ",fontsize=14)
plt.ylabel("BR($N \\rightarrow ll\\nu$) - any flavour neutrino",fontsize=14)
plt.title("Muon-only mixing \n ($x_e$ = 0 , $x_\\mu$ = 1 , $x_\\tau$ = 0)",fontsize=16)
plt.legend(numpoints=1,loc='upper right')

plt.savefig('br_mu_only.pdf')

plt.figure(2)
plt.plot(mass_float,br_IH_eev, marker= 'o',label ="BR($N \\rightarrow ee\\nu$)",markersize=6, color = blue)
plt.plot(mass_float,br_IH_emv, marker= 'v',label ="BR($N \\rightarrow e\\mu\\nu$ )",markersize=6, color = orange)
plt.plot(mass_float,br_IH_mev, marker= 's',label ="BR($N \\rightarrow \\mu e\\nu$)",markersize=4, color = pink)
plt.plot(mass_float,br_IH_mmv, marker= '^',label ="BR($N \\rightarrow \\mu\\mu\\nu$ )",markersize=6, color = purple)
plt.xlabel("$m_N$ [GeV] ",fontsize=14)
plt.ylabel("BR($N \\rightarrow ll\\nu$) - any flavour neutrino",fontsize=14)
plt.title("Inverted hierarchy mixing \n ($x_e$ = 1/3 , $x_\\mu$ = 1/3 , $x_\\tau$ = 1/3)",fontsize=16)
plt.legend(numpoints=1,loc='upper right')
plt.savefig('br_IH.pdf')

plt.figure(3)
plt.plot(mass_float,br_NH_eev, marker= 'o',label ="BR($N \\rightarrow ee\\nu$)",markersize=6, color = blue)
plt.plot(mass_float,br_NH_emv, marker= 'v',label ="BR($N \\rightarrow e\\mu\\nu$ )",markersize=6, color = orange)
plt.plot(mass_float,br_NH_mev, marker= 's',label ="BR($N \\rightarrow \\mu e\\nu$)",markersize=6, color = pink)
plt.plot(mass_float,br_NH_mmv, marker= '^',label ="BR($N \\rightarrow \\mu\\mu\\nu$ )",markersize=6, color = purple)
plt.xlabel("$m_N$ [GeV] ",fontsize=14)
plt.ylabel("BR($N \\rightarrow ll\\nu$) - any flavour neutrino",fontsize=14)
plt.title("Normal hierarchy mixing \n ($x_e$ = 0.06 , $x_\\mu$ = 0.48 , $x_\\tau$ = 0.46)",fontsize=16)
plt.legend(numpoints=1,loc='upper right')
plt.savefig('br_NH.pdf')