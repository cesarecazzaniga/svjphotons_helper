#################################################################################################
## Package to generate self-consistent configuration cards for the pythia 8 hidden valley module for SVJGamma signatures.
## Based on arXiv:2407.08276. Please cite this paper when using this package, in addition to the theory
## papers included in the "citation" string and the software reference.
##
## Written by Cesare Cazzaniga (cesare.cazzaniga@cern.ch)
## adapted from: https://github.com/cms-svj/SVJProduction , https://gitlab.com/simonknapen/dark_showers_tool
#################################################################################################


import os, math, sys, shutil
from string import Template
from glob import glob
import numpy as np
from datetime import date
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

class quark(object):
    def __init__(self,id,mass,charge):
        self.id = id
        self.mass = mass
        self.charge = charge
        self.massrun = mass
        self.color = 3
        self.bf = 1           #theory branching ratios
        self.bf_scaled = 1    #scaled branching ratios by rinv value
        self.on = True        #phase space allowed decay
        self.active = True    #for running nf

    def __repr__(self):
        return str(self.id)+": m = "+str(self.mass)+", mr = "+str(self.massrun)+", on = "+str(self.on)+", bf = "+str(self.bf)

#Added leptons class

class lepton(object):
    def __init__(self,id,mass,charge):
        self.id = id
        self.mass = mass
        self.charge = charge
        self.color = 1
        self.bf = 1           #theory branching ratios
        self.bf_scaled = 1    #scaled branching ratios by rinv value
        self.on = True        #phase space allowed decay  
  
    def __repr__(self):
        return str(self.id)+": m = "+str(self.mass)+", on = "+str(self.on) +", bf = "+str(self.bf)


# follows Ellis, Stirling, Webber calculations
class massRunner(object):
    def __init__(self):
        # QCD scale in GeV
        self.Lambda = 0.218

    # RG terms, assuming nc = 3 (QCD)
    def c(self): return 1./math.pi
    def cp(self,nf): return (303.-10.*nf)/(72.*math.pi)
    def b(self,nf): return (33.-2.*nf)/(12.*math.pi)
    def bp(self,nf): return (153.-19.*nf)/(2.*math.pi*(33.-2.*nf))
    def alphaS(self,Q,nf): return 1./(self.b(nf)*math.log(Q**2/self.Lambda**2))

    # derived terms
    def cb(self,nf): return 12./(33.-2.*nf)
    def one_c_cp_bp_b(self,nf): return 1.+self.cb(nf)*(self.cp(nf)-self.bp(nf))

    # constant of normalization
    def mhat(self,mq,nfq):
        return mq/math.pow(self.alphaS(mq,nfq),self.cb(nfq))/self.one_c_cp_bp_b(nfq)

    # mass formula
    def m(self,mq,nfq,Q,nf):
        # temporary hack: exclude quarks w/ mq < Lambda
        alphaq = self.alphaS(mq,nfq)
        if alphaq < 0: return 0
        else: return self.mhat(mq,nfq)*math.pow(self.alphaS(Q,nf),self.cb(nf))*self.one_c_cp_bp_b(nf)

    # operation
    def run(self,quark,nfq,scale,nf):
        # run to specified scale and nf
        return self.m(quark.mass,nfq,scale,nf)

class quarklist(object):
    def __init__(self):
        # mass-ordered
        self.qlist = [
            quark(2,0.0023,0.67), # up
            quark(1,0.0048,0.33), # down
            quark(3,0.095,0.33),  # strange
            quark(4,1.275,0.67),  # charm
            quark(5,4.18,0.33),   # bottom
        ]
        self.scale = None
        self.runner = massRunner()

    def set(self,scale):
        self.scale = scale
        # mask quarks above scale
        for q in self.qlist:
            # for decays
            if scale is None or 2*q.mass < scale: q.on = True
            else: q.on = False
            # for nf running
            if scale is None or q.mass < scale: q.active = True
            else: q.active = False
        # compute running masses
        if scale is not None:
            qtmp = self.get(active=True)
            nf = len(qtmp)
            for iq,q in enumerate(qtmp):
                q.massrun = self.runner.run(q,iq,scale,nf)
        # or undo running
        else:
            for q in self.qlist:
                q.massrun = q.mass

    def reset(self):
        self.set(None)

    def get(self,active=False):
        return [q for q in self.qlist if (q.active if active else q.on)]



class svjHelper(object):
    def __init__(self,card_author,nevents):
       
        self.quarks_pseudo = quarklist()
        self.quarks_vector = quarklist()
        self.quarks = quarklist()

        # parameters of Hidden sector
        self.n_c = 3
        self.n_f = 2
        
        #card author
        self.username = card_author
        
        #n events to generate
        self.events = nevents


    # Calculation of rho mass from Lattice QCD fits (arXiv:2203.09503v2)
    def calcLatticePrediction(self,mPiOverLambda,mPseudo):        
        mVectorOvermPseudo = (1.0/mPiOverLambda)*math.pow(5.76 + 1.5*math.pow(mPiOverLambda,2) ,0.5) 
        mVector = mVectorOvermPseudo*mPseudo
        return mVector


    # has to be "lambdaHV" because "lambda" is a keyword
    def setModel(self,mMediator,rinv,mPiOverLambda,lambdaHV,brGamma):


        # store the basic parameters
        self.mMediator = mMediator
        self.mPiOverLambda = mPiOverLambda
        self.lambdaHV = lambdaHV
        self.mPseudo = self.mPiOverLambda * self.lambdaHV
        self.mVector = self.calcLatticePrediction(self.mPiOverLambda,self.mPseudo)
        self.rinv = rinv
        self.brGamma = brGamma
    
        # get more parameters
        self.mMin = self.mMediator-1
        self.mMax = self.mMediator+1
        self.mSqua = self.lambdaHV + 0.2 # dark quark constituent mass

        # get limited set of quarks for decays (check mDark against quark masses, compute running)
        self.quarks_pseudo.set(self.mPseudo)
        self.quarks_vector.set(self.mVector)
            
                
    #naming of the output datacard
    def getOutName(self):

        _outname = "SVJGamma"
        _outname += "_mZprime-{:g}".format(self.mMediator)
        _outname += "_rinv-{:g}".format(self.rinv)
        _outname += "_mPioverLambda-{:g}".format(self.mPiOverLambda)
        _outname += "_lambdaHV-{:g}".format(self.lambdaHV)
        _outname += "_brGamma-{:g}".format(self.brGamma)
            
        return _outname

    #implementation of rinv
    def invisibleDecay(self,mesonID,dmID):
        lines = ['{:d}:oneChannel = 1 {:g} 0 {:d} -{:d}'.format(mesonID,self.rinv,dmID,dmID)]
        return lines

    #Here mass insertion for Z' mediated decays to quarks, while br_gamma is used to set the branching ratio to photons for ALP mediated decays 
    def pseudo_scalar_visibleDecay(self,type,mesonID):
        
        theQuarks = self.quarks_pseudo.get()

        if type=="ALPplusMassInsertion":
            denom = sum([q.massrun**2 for q in theQuarks])
            for q in theQuarks:
                q.bf = (1-self.brGamma)*(1.0-self.rinv)*(q.massrun**2)/denom
        else:
            raise ValueError("unknown visible decay type: "+type)
    

        if (type=="ALPplusMassInsertion"):
            lines = ['{:d}:addChannel = 1 {:g} 91 {:d} -{:d}'.format(mesonID,q.bf,q.id,q.id) for q in theQuarks if q.bf>0]
            lines += ['{:d}:addChannel = 1 {:g} 91 22 22'.format(mesonID,self.brGamma*(1.0-self.rinv))]
            

        return lines



   
    ### Vector visible decays via Z' mediated decays
    def vector_visibleDecay(self,type,mesonID):
            theQuarks = self.quarks_vector.get()
           
            if type=="Democratic":
                bfQuarks = (1.0-self.rinv)/float(len(theQuarks))
                for iq,q in enumerate(theQuarks):
                    theQuarks[iq].bf = bfQuarks
            else:
                raise ValueError("unknown visible decay type: "+type)
            lines = ['{:d}:addChannel = 1 {:g} 91 {:d} -{:d}'.format(mesonID,q.bf,q.id,q.id) for q in theQuarks if q.bf>0]
            return lines


    ### Internal decays rho to pi pi
    def vector_internalDecay(self,type,mesonID):
        br = 0 
        if type=="Internal_simple":
            if (mesonID == 4900113):
               br = 1.0  
               decay_prod_1 = 4900211
               decay_prod_2 = -4900211
            if (mesonID == 4900213):
               br = 1.0  
               decay_prod_1 = 4900111
               decay_prod_2 = 4900211
        # lines for decays to quarks                                                                                                                              
        lines_rho_to_pipi = [
             '{:d}:mayDecay=on'.format(mesonID),
             '{:d}:oneChannel = 1 {:g} 91 {:d} -{:d}'.format(mesonID,br,decay_prod_1,decay_prod_2),
        ]
        
        # lines for decays to leptons                                                                                                                              
        return lines_rho_to_pipi



    def getPythiaSettings(self):

        today = date.today()

        #Intro lines
        lines_intro = [
        """##############################################################################################""",
        """# Pythia 8 card for SVJGamma model, as defined in arXiv:2407.08276""",
        f"""# Z mass scale: {self.mMediator} GeV""",
        f"""# rinv: {self.rinv}""",
        f"""# Number colors: {self.n_c}""",
        f"""# Number flavors: {self.n_f}""",
        f"""# Scalar meson mass: {self.mPseudo} GeV""",
        f"""# Vector meson mass: {self.mVector} GeV""",
        f"""# Confinement scale: {self.lambdaHV} GeV""",
        f"""# Branching ratio to photons: {self.brGamma}""",
        f"""# Generated by {' '.join(self.username)} on {today.strftime("%B %d, %Y")}""",
        """# lines starting with ! or # are commented out""",
        """##############################################################################################""",
        '',
        '',
        ]


        #general settings
        lines_settings = [
        '! 1) Settings used in the main program.',
        'Main:numberOfEvents = {}'.format(self.events) ,
        'Main:timesAllowErrors = 3',
        'Random:setSeed = on',
        'Random:seed = 0',
        'Init:showChangedSettings = on',
        'Init:showAllSettings = on',
        'Init:showChangedParticleData = on',
        'Init:showAllParticleData = on',
        'Next:numberCount = 1000',
        'Next:numberShowLHA = 1',
        'Next:numberShowInfo = 1',
        'Next:numberShowProcess = 1',
        'Next:numberShowEvent = 1',
        'Stat:showPartonLevel = on',
        '',
        '',
        ]
        
        #lines beams settings
        lines_beams_settings = [
            '! 2) Beam parameter settings. Values below agree with default ones.',
            'Beams:idA = 2212',
            'Beams:idB = 2212',
            'Beams:eCM = 13600.'
            '',
            '',
        ]

        #Lines s-channel
        lines_schan = [
            '! 3) Turn on the production process.',
            #HV process
            'HiddenValley:ffbar2Zv = on',
            # parameters for leptophobic Z'
            '4900023:m0 = {:g}'.format(self.mMediator),
            '4900023:mMin = {:g}'.format(self.mMin),
            '4900023:mMax = {:g}'.format(self.mMax),
            '4900023:mWidth = 0.01',
            '4900023:doForceWidth = on',
            '',
            '4900023:oneChannel = 1 0.982 102 4900101 -4900101',
            # SM quark couplings needed to produce Zprime from pp initial state
            '4900023:addChannel = 1 0.003 102 1 -1',
            '4900023:addChannel = 1 0.003 102 2 -2',
            '4900023:addChannel = 1 0.003 102 3 -3',
            '4900023:addChannel = 1 0.003 102 4 -4',
            '4900023:addChannel = 1 0.003 102 5 -5',
            '4900023:addChannel = 1 0.003 102 6 -6',
            '',
            '! 4) For MC efficiency let zprime decay only to dark quarks.',
            '4900023:onMode = off',
            '4900023:onIfAny = 4900101',
            '',
            '',
        ]


        #Setting hidden spectrum
        # hidden spectrum:
        # fermionic dark quark,
        # diagonal pseudoscalar meson, off-diagonal pseudoscalar meson, DM stand-in particle,
        # diagonal vector meson, off-diagonal vector meson, DM stand-in particle


        # HV params
        lines_hv = [
            '! 5) HiddenValley Settings',
            'HiddenValley:Ngauge = {:d}'.format(self.n_c),
            'HiddenValley:nFlav = {:d}'.format(self.n_f),
            # when Fv has spin 0, qv spin fixed at 1/2
            'HiddenValley:spinFv = 0',
            'HiddenValley:FSR = on',
            'HiddenValley:fragment = on',
            'HiddenValley:alphaOrder = 1',
            'HiddenValley:setLambda = on',
            'HiddenValley:Lambda = {:g}'.format(self.lambdaHV),
            'HiddenValley:pTminFSR = {:g}'.format(1.1*self.lambdaHV),
            'HiddenValley:probVector = 0.75',
            '',
            '',
        ]
        
        lines_hv.extend([
                '! 6) Dark sector mass spectrum',
                '4900101:m0 = {:g}'.format(self.mSqua),
                '4900111:m0 = {:g}'.format(self.mPseudo),
                '4900211:m0 = {:g}'.format(self.mPseudo),
                '51:m0 = 0.0',
                '51:isResonance = false',
                '4900113:m0 = {:g}'.format(self.mVector),
                '4900213:m0 = {:g}'.format(self.mVector),
                '53:m0 = 0.0',
                '53:isResonance = false',
                # decouple
                '4900001:m0 = 50000',
                '4900002:m0 = 50000',
                '4900003:m0 = 50000',
                '4900004:m0 = 50000',
                '4900005:m0 = 50000',
                '4900006:m0 = 50000',
                '4900011:m0 = 50000',
                '4900012:m0 = 50000',
                '4900013:m0 = 50000',
                '4900014:m0 = 50000',
                '4900015:m0 = 50000',
                '4900016:m0 = 50000',
                '',
                '',
                '! 7) HV decays',
        ])


        # branching - effective rinv (applies to all meson species)
        # pseudoscalars have mass insertion decays (from Z') and ALP decays
        # vectors have democratic decays from Z' or internal decays to pi pi
      
        #all pseudoscalars decays
        lines_hv += self.invisibleDecay(4900111,51)
        lines_hv += self.pseudo_scalar_visibleDecay("ALPplusMassInsertion",4900111)
        lines_hv += self.invisibleDecay(4900211,51)
        lines_hv += self.pseudo_scalar_visibleDecay("ALPplusMassInsertion",4900211)
        
        if self.mPiOverLambda <= 1.5:
            #all vector mesons decays
            lines_hv += self.vector_internalDecay("Internal_simple",4900113)
            lines_hv += self.vector_internalDecay("Internal_simple",4900213)
        else:
            lines_hv += self.invisibleDecay(4900113,53)
            lines_hv += self.vector_visibleDecay("Democratic",4900113)
            lines_hv += self.invisibleDecay(4900213,53)
            lines_hv += self.vector_visibleDecay("Democratic",4900213)


        lines = lines_intro + lines_settings + lines_beams_settings + lines_schan + lines_hv

        return lines

    def get_citations(self):
         
        citation_string = """
         @software{Cazzaniga_SVJGamma_models_production_2024,
            author = {Cazzaniga, Cesare},
            month = dec,
            title = {{SVJGamma models production}},
            url = {https://github.com/cesarecazzaniga/svjphotons_helper},
            version = {1},
            year = {2024}
         }
         
         """
         
            
        citation_string+="""@article{Cazzaniga:2024mmv,
        author = \"Cazzaniga, Cesare and Russo, Alessandro and Sitti, Emre and de Cosa, Annapaola\",
        title = \"{Phenomenology of photons-enriched semi-visible jets}\",
        eprint = \"2407.08276\",
        archivePrefix = \"arXiv\",
        primaryClass = \"hep-ph\",
        doi = \"10.1140/epjc/s10052-024-13613-9\",
        journal = \"Eur. Phys. J. C\",
        volume = \"84\",
        number = \"11\",
        pages = \"1223\",
        year = \"2024\"
    }
        
        @article{Knapen:2021eip,
        author         = \"Knapen, Simon and Shelton, Jessie and Xu\",
        title          = \"{Perturbative benchmark models for a dark shower search program}\",
        journal        = \"Phys. Rev. D\",
        volume         = \"103\",
        year           = \"2021\",
        pages          = \"115013\",
        doi            = \"10.1103/PhysRevD.103.115013\",
        eprint         = \"2103.01238\",
        archivePrefix  = \"arXiv\",
        primaryClass   = \"hep-ph\"
    }
        
        """
        
            
        return citation_string

    
if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("--mZprime", dest="mZprime",  required=True, type=float, help="Zprime mass (GeV)")
    parser.add_argument("--rinv", dest="rinv", required=True, type=float, help="invisible fraction")
    parser.add_argument("--mPiOverLambda", dest="mPiOverLambda", required=True, type=float, help="Lightest dark pseudoscalar mass over LambdaHV")
    parser.add_argument("--lambda", dest="lambdaHV", required=True, type=float, help="dark sector confinement scale")
    parser.add_argument("--brGamma", dest="brGamma", required=True, type=float, help="Branching to photons of dark pions")
    parser.add_argument("--nevents", dest="nevents", required=True, type=int, help="Number of events to generate")
    parser.add_argument("--card_author", dest="card_author", required=False, type=str,nargs='+', default="", help="author of the generated card")
    args = parser.parse_args()


    helper = svjHelper(args.card_author, args.nevents)
    helper.setModel(args.mZprime, args.rinv, args.mPiOverLambda, args.lambdaHV, args.brGamma)

    lines = helper.getPythiaSettings()
    fname = helper.getOutName()+".cmnd"
    with open(fname,'w') as file:
        file.write('\n'.join(lines))
        
    print("Generated card and saved in file: ", fname)
    
    #print citation string
    citations = helper.get_citations()
    print("Please cite: ")
    print(citations)
    
        
        
        
