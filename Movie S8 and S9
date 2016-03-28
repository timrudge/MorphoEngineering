import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator #add
import numpy
import math

#Specify parameter for solving diffusion dynamics 
grid_dim = (64, 8, 12) # dimension of diffusion space, unit = number of grid
grid_size = (4, 4, 4) # grid size
grid_orig = (-128, -14, -8) # where to place the diffusion space onto simulation space

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False, max_cells=30000,max_sqs=192**2)

    regul = ModuleRegulator(sim, sim.moduleName)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    sim.addCell(cellType=0, pos=(0,0,0))

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = 10

def numSpecies():
    return(0)

def numSignals():
    return(0)

def init(cell):
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    cell.growthRate = 1.0
    cell.n_a = 1 
    cell.n_b = 1

def update(cells):
    for (id, cell) in cells.iteritems():
        if len(cells)>1000:
            gr1 = 1.0
            gr2 = 0.5
            cell.growthRate = (gr2-gr1)*cell.n_a*0.5 + gr1

        cell.color = [0.1, cell.n_a/2.0, cell.n_b/2.0]
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    d1.targetVol = 3.5 + random.uniform(0.0,0.5)
    d2.targetVol = 3.5 + random.uniform(0.0,0.5)
    plasmids = [0]*parent.n_a*2 + [1]*parent.n_b*2
    random.shuffle(plasmids)
    d1.n_a = 0
    d1.n_b = 0
    d2.n_a = 0
    d2.n_b = 0
    for p in plasmids[:2]:
        if p == 0: d1.n_a +=1
        else: d1.n_b +=1
    for p in plasmids[2:4]:
        if p == 0: d2.n_a +=1
        else: d2.n_b +=1
    assert parent.n_a + parent.n_b == 2
    assert d1.n_a + d1.n_b == 2
    assert d2.n_a + d2.n_b == 2
    assert parent.n_a*2 == d1.n_a+d2.n_a
    assert parent.n_b*2 == d1.n_b+d2.n_b
    assert parent.n_a > 0 or (d1.n_a == 0 and d2.n_a == 0)
    assert parent.n_b > 0 or (d1.n_b == 0 and d2.n_b == 0)
