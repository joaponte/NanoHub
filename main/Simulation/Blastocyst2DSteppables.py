from PySteppables import *
#from Parameters import *
import CompuCell
import sys
from random import randint
import CompuCellSetup
from PySteppablesExamples import MitosisSteppableBase
from PlayerPython import *
from math import *
import numpy as np

#Time for division
TCycle=1400 #1200
TRelax = 200
scale = 1.0
scale_update = 0.9

#Fluid
rr=15 #Radius for random location Sacled to lattice dimensions
fractionfluidmax=.25 #Max fraction of interior area occupied by fluid
pressure_factor = 0.0 #pressure factor

#GRN
sigma = 0.45
iteration = 1
update_frequency = 10.0
dt = 0.01 * update_frequency
h = 4
import global_vars

class Original_Config(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        global_vars.r1= self.dim.x/2+1#ZP External radius
        global_vars.r2 = 0.9*(self.dim.x/2+1)#ZP Internal Radius
        global_vars.r3 = global_vars.r2/1.125 #Initial Cell Radius #1.2
        
        cella=self.newCell(self.PELLUCIDA)
        cellb=self.newCell(self.BLASTOMERE)
        for x in range(0,self.dim.x):
            for y in range(0,self.dim.y):
                rd=sqrt((x-(self.dim.x/2))**2+(y-(self.dim.y/2))**2)
                if rd<=global_vars.r3:
                    self.cellField[x,y,0]=cellb
                elif(global_vars.r2<rd<global_vars.r1):
                    self.cellField[x,y,0]=cella
        
        for cell in self.cellListByType(self.BLASTOMERE):
            cell.targetVolume= pi*(global_vars.r3)**2 #Setting tgVol to be actual volume
            cell.lambdaVolume=1.0 #Inverse compressibility (is it known? how its related? effect of change
            cell.dict['Time'] = 600
            cell.dict['Boundary'] = 1


class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
    
    def start(self):
        self.time_after_division = 0
        
    def step(self,mcs):
        for cell in self.cellListByType(self.BLASTOMERE,self.TROPHOBLAST,self.INNERCELL):
            cell.dict['Time'] -= 1
            if cell.dict['Time'] == 0:
                if global_vars.gen == 1:
                    self.divideCellOrientationVectorBased(cell, 1, 0, 0)
                elif global_vars.gen == 2:
                    xc = cell.xCOM-(self.dim.x/2)
                    yc = cell.yCOM-(self.dim.y/2)
                    L=sqrt(xc**2+yc**2)
                    self.divideCellOrientationVectorBased(cell, yc/L, - xc/L, 0)
                else:
                    self.divideCellRandomOrientation(cell)
        
        num_cells = len(self.cellListByType(self.BLASTOMERE,self.TROPHOBLAST,self.INNERCELL))
        if  num_cells == 2 ** (global_vars.gen) and global_vars.gen < 6:
            self.time_after_division += 1
            if self.time_after_division == TRelax:
                global_vars.gen += 1
                for cell in self.cellListByType(self.BLASTOMERE,self.TROPHOBLAST,self.INNERCELL):
                    cell.dict['Time'] = 50*randint(TCycle/50, (TCycle + 2*TRelax)/50)
        else:
            self.time_after_division = 0
    
    def updateAttributes(self):
        self.parentCell.targetVolume /= 2.0 # reducing parent target volume by half
        self.cloneParent2Child()

class Fluid_Growth(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        self.fluid_growth_rate = 0
        self.flag2 = 0
        
    def step(self,mcs):
        if global_vars.gen ==  5 and self.flag2 ==0:
            flag1 = 0
            while flag1 == 0:
                #x=randint((self.dim.x/2)-rr,(self.dim.x/2)+rr)
                x=randint(60,70)
                #y=randint((self.dim.y/2)-rr,(self.dim.y/2)+rr)
                y=randint(60,70)
                cell = self.cellField[x, y, 0]
                if cell:
                    pass
                else:
                    flag1 = 1
            maxfluidsize = fractionfluidmax*pi*global_vars.r2**2
            inifluidsize = 0.025*maxfluidsize
            self.fluid_growth_rate = (maxfluidsize-inifluidsize)/(12000-mcs)
            Fluid=self.newCell(self.FLUID)
            self.cellField[x:x+sqrt(inifluidsize),y:y+sqrt(inifluidsize),0]=Fluid
            Fluid.targetVolume = inifluidsize
            Fluid.lambdaVolume = 1.0
            self.flag2 = 1
            
        ncell = len(self.cellListByType(self.BLASTOMERE,self.TROPHOBLAST,self.INNERCELL))
        if global_vars.gen >= 5:
            for cell in self.cellListByType(self.FLUID):
                cell.targetVolume+=self.fluid_growth_rate
                if global_vars.gen == 6:
                    cell.lambdaVolume = 0.1
            for cell in self.cellListByType(self.BLASTOMERE,self.TROPHOBLAST,self.INNERCELL):
                cell.targetVolume-=self.fluid_growth_rate/ncell * pressure_factor
        

class Boundary(SteppableBasePy):        
    def start(self):
        self.Zeta = 0.5
        self.rb = 42
        
    def step(self,mcs):
        self.dAdhesion = 2.5*self.Zeta   
        if global_vars.gen >= 5 and mcs%update_frequency == 0: #5
            for cell in self.cellListByType(self.BLASTOMERE,self.TROPHOBLAST,self.INNERCELL):
                cell.dict['Boundary'] = 1
                xCOM = cell.xCOM
                yCOM = cell.yCOM
                dist_from_center = sqrt((xCOM-self.dim.x/2)**2 + (yCOM-self.dim.y/2)**2)
                total_membrane_pix = 0.0
                no_contact_pix = 0.0
                if dist_from_center >= global_vars.r2 - global_vars.r3:
                    pixelList = self.getCellBoundaryPixelList(cell)
                    for pixel in pixelList:
                        total_membrane_pix+=1
                        pix_dist_from_center=sqrt((pixel.pixel.x-(self.dim.x/2))**2+(pixel.pixel.y-(self.dim.y/2))**2)
                        #rb = gloabl_vars.r2 - self.Zeta
                        if self.rb>=pix_dist_from_center:
                            no_contact_pix+=1
                    cell.dict['Boundary'] = no_contact_pix/total_membrane_pix
                    #print cell.dict['Boundary']


class GRN(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        for cell in self.cellListByType(self.BLASTOMERE):
            cell.dict["Cdx2"] = 0.0
            cell.dict["Oct4"] = 0.0      

    def step(self,mcs):
        def Non_Dimensional_GRN(x,y,cell):
            k = 0.7
            b = 0.7
            I = 0.6
            S = 1.3 - 1.5 * cell.dict['Boundary']
            dx = k * ( b + S + x**h / ( x**h + 0.5**h ) ) * ( ( 1 - I ) + I * 0.5**h / ( y**h + 0.5**h ) ) - x
            dy = k * ( b + y**h / ( y**h + 0.5**h ) ) * ( ( 1 - I ) + I * 0.5**h / ( x**h + 0.5**h ) ) - y
            return [dx,dy]
        
        def Dimensional_GRN(x,y,cell):
            ax=1
            ay = 1 
            delta=1
            Ix=.4
            Iy=.4
            by=.7
            k = 0.8
            bx = 2-1.5*cell.dict['Boundary']
            dx=k*(bx+ax*x**h/(0.5**h+x**h))*(Ix+(1-Ix)*0.5**h/(0.5**h+y**h))-x;
            dy=k*(by+ay*y**h/(0.5**h+y**h))*(Iy+(1-Iy)*0.5**h/(0.5**h+x**h))-delta*y
            return [dx,dy]
            
        if global_vars.gen >= 5 and mcs%update_frequency == 0:
            for cell in self.cellListByType(self.BLASTOMERE,self.TROPHOBLAST,self.INNERCELL):
                #Get values at previous step
                x = cell.dict["Cdx2"]
                y = cell.dict["Oct4"]
                #ODES
                #ODE = Non_Dimensional_GRN(x,y,cell)
                ODE = Dimensional_GRN(x,y,cell)
                dx = ODE[0]
                dy = ODE[1]
                #Step forward
                x1 = x + dt * dx + sqrt(dt) * np.random.normal() * sigma * x
                y1 = y + dt * dy + sqrt(dt) * np.random.normal() * sigma * y
                #Determines cell type
                if  0.8 < x1 - y1:
                    cell.type = self.TROPHOBLAST
                elif 0.8 < y1 - x1:
                    cell.type = self.INNERCELL
                #Update values
                cell.dict["Cdx2"] = x1
                cell.dict["Oct4"] = y1

class GRN_Field(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.Cdx2Field=CompuCellSetup.createScalarFieldCellLevelPy("Cdx2")
        self.Oct4Field=CompuCellSetup.createScalarFieldCellLevelPy("Oct4")
    def step(self,mcs):
        self.Cdx2Field.clear()
        self.Oct4Field.clear()
        for cell in self.cellListByType(self.BLASTOMERE,self.TROPHOBLAST,self.INNERCELL):
            self.Cdx2Field[cell]=cell.dict["Cdx2"]
            self.Oct4Field[cell]=cell.dict["Oct4"]
