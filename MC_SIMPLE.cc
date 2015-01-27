// -*- C++ -*-
//Mundane include stuff, no need to bother about usually

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
	
//Necessary declarations, again no need to play with 

namespace Rivet {

  using namespace Cuts;

  // A very simple analysis looking at some final state objects
  // Charged particles, jets, Z bosons
  // Again, just need to change the analysis name if needed

  class MC_SIMPLE : public Analysis {
  public:

    /// Constructor
    MC_SIMPLE()
      : Analysis("MC_SIMPLE")
    {    }


  public:

 

    void init() {

    // Declare the projections, needed one for accessing each type object we will use
    // Lets look at charged particles and jets ;-)
    // Note that cfs, fs etc are just named...

    // Charged final state for the distributions, within |eta| < 2.5 and pT > 0.5 GeV. 
     const ChargedFinalState cfs(-2.5, 2.5, 0.5*GeV);
     addProjection(cfs, "CFS");

    // For jets, we need inputs for jet algorithm. We choose all particles within |eta| < 5.0
      const FinalState fs(-5.0, 5.0, 0.0*GeV);
      addProjection(fs, "FS");
    //Fastjet is the separate code which makes jets, here using ANTIKT algorithm with radius 0.4
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets");

     //We can add more projections if needed, lets say for Z or W boson, more types of jets, more types of particles etc
 
    // Declare the histograms. Add more as needed, keeping in mind the type

      _hChargedEta = bookHisto1D("Eta",50,-5,5);
      _hSumpT = bookHisto1D("Sum pT",100,0,100);
      _hNJets = bookHisto1D("N Jets",10,0,10);
      _hJetpT = bookHisto1D("Jet pT",100,0,100);

    } // end of init


  // First loop is over all events, one by one. We also save event weight as a good practise.

    void analyze(const Event& event) {
      const double weight = event.weight();

      //Start by looking at charged particles, apply the previously declared projection
      const FinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
    
      //Two methods to access particles, first by saving them in vector "particles", where they can be accesed by indices.
      Particles particles = cfs.particlesByPt();
      if(particles.size()>0) {
        const double pTlead = particles[0].pT();
        cout << " Just to test, pT of leading particle: "<< pTlead << endl;
      }

      double pTsum =0; //Declare and initialize
      // Loop over all (charged) particles, each indicated by "p". Each p is a particle at the loop poistion.
      foreach (const Particle& p, cfs.particles()) {

	//Read off eta and pT of all p's
        const double pT = p.pT();
        const double eta = p.eta();
	//We want to calculate sum of all particle pT in each event
        pTsum = pTsum + pT;

	//Fill eta histogram inside the particle loop, since its for all particles
        _hChargedEta->fill(eta, weight);

      }//End loop over particles

      //Fill Sum pT after summing over all particles 
      _hSumpT->fill(pTsum, weight);
  
      //Now look at jets, first apply the previously declared projection
      const FastJets& fastjets = applyProjection<FastJets>(event, "Jets");  

      //Method 1 to access the jets, by using array of "jets"
      const Jets jets = fastjets.jetsByPt(20.*GeV); 
      const double NJ = jets.size();
 
      cout << " Number of jets: " << NJ << endl;   
 
      _hNJets->fill(NJ, weight);
      if(NJ>0) cout << " Leading jet pT: " << jets[0].pT() << endl;

      //Method 2 to access jets, by looping over them, requiring minimum pT of 20 GeV
      //The loop means "j" is the each jet object inside the loop, at loop position

        foreach(Jet j, fastjets.jetsByPt(20.*GeV)) {
	   _hJetpT->fill(j.pT(), weight);

        }//End loop over jets

    }//End loop over events



 // Finalize: this would contain normalizing the histograms.
    void finalize() {}


 private:

    // Save pointers to the Histograms, again separately for each types
    Histo1DPtr _hChargedEta, _hSumpT, _hNJets, _hJetpT;
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_SIMPLE);

}//Grand end
