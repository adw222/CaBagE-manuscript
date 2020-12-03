# CaBagE Protocol
**Samples**: 5 ug high molecular weight genomic DNA

**Guides:** Assembled crRNA+tracr gRNA pair flanking target region in PAM-in orientation

**NEBNext FFPE Repair**
Mix the following components in a sterile nuclease-free tube:
_Reagent						X1_
DNA (200ng/uL)					14 µl (5.01 ug)
IDTE 1X							39.5 µl
FFPE DNA Repair Buffer			6.5 µl
NEBNext FFPE DNA repair Mix		2 µl

Total Volume						62 µl

- [ ] Mix by flicking followed by a quick spin to collect all liquid from the sides of the tube.
- [ ] Incubate at 20°C for 15 minutes.
Cleanup Using AMPure XP Beads
- [ ] Vortex AMPure XP Beads to resuspend.
- [ ] Add 1.8X [111.6uL] of resuspended AMPure XP Beads to the repair reaction. Mix by flicking.
- [ ] Incubate on hula mixer for 5 minutes at room temperature.
- [ ] Put the tube/PCR plate on an appropriate magnetic stand to separate beads from supernatant. After the solution is clear (about 5 minutes), carefully remove and discard the supernatant. Be careful not to disturb the beads that contain the DNA targets.
- [ ] Add 200 μl of 70% freshly prepared ethanol to the tube/PCR plate while in the magnetic stand. Incubate at room temperature for 30 seconds, and then carefully remove and discard 168uL of supernatant.
- [ ] Repeat previous step once.
- [ ] Air dry beads for up to 5 minutes while the tube/PCR plate is on the magnetic stand with the lid open.
Caution: Do not over-dry the beads. This may result in lower recovery of DNA target. Elute the samples when the beads are still dark brown and glossy looking, but when all
visible liquid has evaporated. When the beads turn lighter brown and start to crack, they are too dry.
- [ ] Remove the tube_plate from the magnet. Elute DNA target by adding 20 μl H2O to the beads. Mix well on a vortex mixer or by pipetting up and down, and incubate for 10 minutes at room temperature. Put the tube_PCR plate in the magnetic stand until the solution is clear.
- [ ] Without disturbing the bead pellet, carefully transfer 20 μl of the supernatant to a fresh, sterile microfuge tube.


**Prepare working mix of HiFi Cas9 ice:**
	Dilute the 10X CutSmart buffer 1:8 with Nuclease Free Water (NFW) to 	make 1.25X CutSmart Buffer (14uL NFW + 2uL 10X CutSmart Buffer)

	Dilute the HiFi Cas9 1:5 using the 1.25X CutSmart made above (1uL  	
	HiFi Cas9 + 4uL 1.25X CutSmart)

**Cas9/Guide RNA Pre-Assembly**
Stock						Volume (uL)		Final Concentration
dH2O						19.3			
10X Cutsmart Buffer			2.5				1X
10 uM sgRNA For				0.6				150nM
10 uM sgRNA Rev			0.6				150nM
1:5 diluted 62μM HiFI Cas9	0.5				155nM

Total Volume:				23.5

- [ ] Incubate 10 minutes at 25°C

**Cas9 Digestion**
_Add to Previous Reaction Product [23.5uL]:_
10X Cutsmart Buffer			1.5				1X
15 nM substrate DNA			15				1.5nM

Total Volume:				40

- [ ] Mix thoroughly and pulse-spin in a microfuge
- [ ] Incubate at 37°C for 15 minutes

**While Cas9 is incubating:**

- [ ] Dilute the 10X CutSmart buffer with Nuclease Free Water (NFW) to make 3X CutSmart Buffer
_Stock						Volume (uL) X1_
NFW						7
10X CutSmart Buffer			3

- [ ] Make Exo Mix
_Stock						Volume (uL) X1_
3X CutSmart Buffer			4
Lambda Exo					4
Exo I						2
Exo III						2

**Exonuclease Digestion**
Add to Previous Reaction Product [40uL]:
_Stock						Volume (uL) X1_
Diluted Exo Mix				10

Total Volume:				50

- [ ] Incubate at 37°C for 2 hours, 80°C for 20 minutes

**A-tailing**
Add to Previous Reaction Product [50uL]:
_Stock						Volume (uL)_
10mM dATP					1
Taq Polymerase				1

Total Volume:				52

- [ ] Incubate at 72°C for 5 minutes
- [ ] Hold at 12°C

**During this incubation, get out the following components from SQK-LSK109 reagents to thaw:**
LNB (keep at room temp)
AMX (keeping on ice after thaw)

**Adapter ligation and  Ampure Purification**
- [ ]  In a separate tube, make "Ligation mix”: NOTE: Make separately for each sample if generating multiple libraries
_Stock						Volume (uL) X1	Notes_
LNB							25				ligation buffer from 													LSK109 kit, must mix by 												pipetting because LNB is 												very viscous
NFW						13
Quick T4 Ligase				5				NEB (E6057)
AMX						5				sequencing adaptor from 												LSK109 kit, AMX will 													degrade in this solution, 												add right before use

- [ ] Transfer samples to 1.5mL Eppindorf Tube
- [ ] Using pipette, mix ligation mix gently until homogenous (LNB is very viscous)
- [ ] Add half (24uL) of "ligation mix" to DNA*, mix by gentle flicking
- [ ] Add remaining (~24uL) "ligation mix" to DNA [total volume = 100uL]
- [ ] Mix by gentle flicking
* A DNA precipitate made from at this step, adding "ligation mix" in 2-steps helps to reduce this
* Formation of DNA precipitate does not appear to interfere with protocol efficiency
- [ ] Rotate ligation for 10min on hula mixer at Room Temperature

**During incubation, get out SFB and EB from LSK109 kit and Ampure beads, keeping all at room temp**

- [ ] After ligation, add equivolume (100uL) amount of IDTE (pH 7.5) [new final volume = 200uL ]
- [ ] Add 0.3X Ampure (60uL if DNA is currently in 200uL) and mix by gentle flicking until homogeneous
- [ ] Rotate for 5 minutes on a hula mixer, then keep on benchtop for 5 minutes to allow Ampure beads to bind DNA
- [ ] Place on magnetic tube rack - allow 2.5 minutes for beads to collect on back of tube
- [ ] Remove supernatant with a pipette, taking care not to disturb the beads
- [ ] Add 200uL of SFB to beads. Remove tube from magnetic stand and resuspend by gentle flicking
- [ ] Return sample to magnetic tube rack , allowing 2.5 minutes for beads to collect
- [ ] Remove supernatant
- [ ] Repeat steps for a second bead wash w/ 200uL SFB
- [ ] Remove supernatant, briefly spin tube to collect any excess SFB and carefully remove remaining solution from the beads
- [ ] Add 16.6uL of elution buffer (EB), mix by flicking, and keep at room temperature for 10 minutes to elute
- [ ] Return tube with sample to magnetic rack , and allow 2min for beads to collect
- [ ] Transfer (and keep!) ~15.8uL supernatant with DNA to a new PCR tube

**Library Prep and Sequencing**

- [ ] Draw a few uLs fluid form priming port of FLO-MIN106
- [ ] Perform initial priming of flow cell by loading 800uL FLB into MinION priming port - wait 5 minutes
- [ ] Finish library prep by adding the following components to the DNA eluate:
26uL SQB
9.5uL Loading beads (LB)
0.5uL SQT
- [ ] Perform second flow cell priming with 200uL of a [70uL SQB + 70uL NFW + 70uL FLB], then immediately afterwards load sample dropwise into MinION sample loading port
- [ ] Start sequencing run using MinKNOW software
