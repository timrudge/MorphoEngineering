import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
from CellModeller.GUI import Renderers
from CellModeller.Integration.CLCrankNicIntegrator import CLCrankNicIntegrator #add
import numpy
import math


def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, jitter_z=False, max_cells=50000)

    regul = ModuleRegulator(sim, sim.moduleName)	# use this file for reg too
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    sim.addCell(cellType=0, pos=(0,0,0))

    # Add some objects to draw the models
    therenderer = Renderers.GLBacteriumRenderer(sim)
    sim.addRenderer(therenderer)
    sim.pickleSteps = 10

def init(cell):
    cell.targetVol = 3.5 + random.uniform(0.0,0.5)
    cell.growthRate = 1.0
    cell.n_a = 1 
    cell.n_b = 1

def init(cell):

    cell.targetVol = 3.5 + random.uniform(0.0,0.5)

    cell.growthRate = 1.0

    cell.n = [1,1,1] 



def update(cells):

    for (id, cell) in cells.iteritems():

        r = 0.3+cell.n[0]/3.0

        g = 0.3+cell.n[1]/3.0

        b = 0.3+cell.n[2]/3.0

        cell.color = [r,g,b]

        if cell.volume > cell.targetVol:

            cell.divideFlag = True



def divide(parent, d1, d2):

    d1.targetVol = 3.5 + random.uniform(0.0,0.5)

    d2.targetVol = 3.5 + random.uniform(0.0,0.5)

    plasmids = []

    for i in range(3):

        plasmids += [i]*parent.n[i]*2

        d1.n[i] = 0

        d2.n[i] = 0

    random.shuffle(plasmids)

    for p in plasmids[:3]:

        d1.n[p] += 1

    for p in plasmids[3:6]:

        d2.n[p] += 1

    '''

    assert parent.n_a + parent.n_b == 2

    assert d1.n_a + d1.n_b == 2

    assert d2.n_a + d2.n_b == 2

    assert parent.n_a*2 == d1.n_a+d2.n_a

    assert parent.n_b*2 == d1.n_b+d2.n_b

    assert parent.n_a > 0 or (d1.n_a == 0 and d2.n_a == 0)

    assert parent.n_b > 0 or (d1.n_b == 0 and d2.n_b == 0)

    '''

