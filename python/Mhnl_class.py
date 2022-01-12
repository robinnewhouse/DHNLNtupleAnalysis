import ROOT
import numpy as np


class Mhnl:
    """
    The class used to compute the HNL mass (m_HNL) for each event.
    Inputs:
    - pv: ROOT TVector3  object for the location (x,y,z) of the primary vertx (PV). See: https://root.cern.ch/doc/master/classTVector3.html
    - dv: ROOT TVector3 object for the location (x,y,z) of the displaced vertx (DV). See: https://root.cern.ch/doc/master/classTVector3.html
    - displaced_leptons: Python list of ROOT TLorentzVector objects (px, py, pz, E) for two tracks in the DV. See: https://root.cern.ch/doc/master/classTLorentzVector.html
    - prompt_lepton: ROOT TLorentzVector object (px, py, pz, E) for prompt lepton. See: https://root.cern.ch/doc/master/classTLorentzVector.html
    - fix_w_mass: A boolean flag that fixes the W mass to the pole mass in the energy momentum conservation expression. Default: False
    """

    def __init__(self, pv, dv, displaced_leptons, prompt_lepton, fix_w_mass=False):
        MW = 80.379  # W boson pole mass in GeV
        MW2 = MW ** 2
        WGamma = 2.085  # W boson width in GeV
        self.mhnl = -1  # HNL mass (m_HNL)
        self.alt_mhnl = -1  # Alternative solution to the quadratic equation
        self.hnl_pt = -1
        self.hnl_eta = -99
        self.hnl_phi = -99
        self.m_lll = -1
        self.no_solution = -1

        # Assign p0 to be the 3-vector of the prompt lepton
        p0 = ROOT.TVector3(prompt_lepton.Px(), prompt_lepton.Py(), prompt_lepton.Pz())

        # Assign d0 and d1 to be the 3-vectors for the displaced leptons
        d0 = ROOT.TVector3(displaced_leptons[0].Px(), displaced_leptons[0].Py(), displaced_leptons[0].Pz())
        d1 = ROOT.TVector3(displaced_leptons[1].Px(), displaced_leptons[1].Py(), displaced_leptons[1].Pz())

        # Choose z direction to be along decay
        decayV = dv - pv
        z = decayV * (1.0 / decayV.Mag())

        # Visible (2 leptons in DV) system
        dvis = d0 + d1
        mvis2 = 2 * (d0.Mag() * d1.Mag() - d0.Dot(d1))

        # Define plane perpendicular to direction of decay
        x = dvis.Cross(z)
        qperp = x.Mag()
        qperp2 = qperp * qperp

        # x and y unit vectors
        x = x * (1.0 / qperp)
        y = z.Cross(x)

        # Visible momentum in the (x,y,z) coordinates
        qv = ROOT.TVector3(0., qperp, dvis.Dot(z))
        Ev = np.sqrt(mvis2 + qv.Dot(qv))

        # Prompt lepton in new (x,y,z) coordinates
        pp = ROOT.TVector3(p0.Dot(x), p0.Dot(y), p0.Dot(z))
        Ep = np.sqrt(pp.Dot(pp))

        # Terms from conservation of 4 momentum that involve various powers of the neutrino z momentum
        B = (pp[2] + qv[2])
        E = Ep + Ev
        B = B / E

        # Minimum possible W mass
        alpha = qperp * B / np.sqrt(1 - B ** 2)
        qn1 = ROOT.TVector3(0, -qperp, alpha)
        En1 = qn1.Mag()
        qtot1 = pp + qv + qn1
        Etot1 = Ep + Ev + En1
        mWMin2 = Etot1 ** 2 - qtot1.Dot(qtot1)

        # Get median of the cumulative distribution
        cdVal = 0.5 + np.arctan((mWMin2 - MW2) / MW / WGamma) / np.pi
        cdMed = (1 + cdVal) / 2

        # Invert to get mass of median allowed range
        mMed2 = MW2 + MW * WGamma * np.tan(np.pi * (cdMed - 0.5))

        # W mass choice:
        if fix_w_mass:
            MW2fit = MW2  # Use pole mass for W mass
        else:
            MW2fit = mMed2  # Use median W mass in kinematically allowed regions

        A = (MW2fit - mvis2) / 2. - Ep * Ev - qperp2 + pp[2] * qv[2]
        A = A / E

        # Coefficients of the quadratic to solve
        b = 2 * A * B / (B * B - 1)
        c = (A * A - qperp2) / (B * B - 1)
        # Arugment in the square root of the quadric solution
        arg = b * b - 4 * c
        # Protect against imaginary solutions
        if arg > 0:
            rad = np.sqrt(arg)
            self.no_solution = 0  # real solution exists for the quadratic equation
        else:
            rad = -1000
            self.no_solution = 1  # no real solution to the quadratic equation

        # These are the possible z momenta for the neutrino
        sol1 = (-b + rad) / 2
        sol2 = (-b - rad) / 2

        # Make 3-vectors of the two z-momentum solutions for the neutrino
        qn1 = ROOT.TVector3(0, -qperp, sol1)
        En1 = qn1.Mag()
        qn2 = ROOT.TVector3(0, -qperp, sol2)
        En2 = qn2.Mag()

        # Total momentum of the decaying system
        qtot1 = qv + qn1
        qtot2 = qv + qn2

        # Neutrino momentum in original coordinates
        dn1 = y * (-1 * qperp) + z * sol1
        dn2 = y * (-1 * qperp) + z * sol2

        # Make 4-vectors in the original coordinates
        # Use pion mass assumption for mass of tracks (default for ATLAS tracking)
        pnu2 = ROOT.TLorentzVector(dn2, dn2.Mag())
        pnu1 = ROOT.TLorentzVector(dn1, dn1.Mag())
        pion_mass = 0.139  # pion mass in GeV
        plep0 = ROOT.TLorentzVector()
        ptrk0 = ROOT.TLorentzVector()
        ptrk1 = ROOT.TLorentzVector()
        plep0.SetPxPyPzE(p0.X(), p0.Y(), p0.Z(), np.sqrt(p0.Mag() ** 2 + pion_mass ** 2))
        ptrk0.SetPxPyPzE(d0.X(), d0.Y(), d0.Z(), np.sqrt(d0.Mag() ** 2 + pion_mass ** 2))
        ptrk1.SetPxPyPzE(d1.X(), d1.Y(), d1.Z(), np.sqrt(d1.Mag() ** 2 + pion_mass ** 2))
        pHNL1 = pnu1 + ptrk0 + ptrk1
        pHNL2 = pnu2 + ptrk0 + ptrk1
        plll = plep0 + ptrk0 + ptrk1

        # Set the final attributes of the class
        self.mhnl = pHNL1.M()
        self.hnl_pt = pHNL1.Pt()
        self.hnl_eta = pHNL1.Eta()
        self.hnl_phi = pHNL1.Phi()
        self.m_lll = plll.M()
        self.alt_mhnl = pHNL2.M()


'''Here is an example event to show the functionality of the above class.'''
if __name__ == "__main__":
    print("HNL mass calculation examples", "\n")

    # example 1
    print("Example 1 is from an MC generated 10 GeV HNL")
    dv_1 = ROOT.TVector3(15.512019157409668, 19.17801856994629, -107.96952056884766)
    pv_1 = ROOT.TVector3(-0.5000697374343872, -0.5047531723976135, -40.7725830078125)
    displaced_leptons_1 = [ROOT.TLorentzVector(), ROOT.TLorentzVector()]
    displaced_leptons_1[0].SetPxPyPzE(2.48017372631826, -0.2136493157850394, -4.773501905451388, 666)
    displaced_leptons_1[1].SetPxPyPzE(9.90227566884341, 15.053740580375933, -43.85763404407506, 666)
    prompt_lepton_1 = ROOT.TLorentzVector()
    prompt_lepton_1.SetPxPyPzE(-19.65923058291866, -32.908382397006065, -198.16919648843634, 666)
    example_1 = Mhnl(pv_1, dv_1, displaced_leptons_1, prompt_lepton_1, fix_w_mass=False)

    print("Truth HNL mass:", "10 GeV")
    print("Calculated HNL mass:", example_1.mhnl, "GeV")
