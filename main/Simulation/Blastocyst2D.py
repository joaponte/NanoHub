
import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])


import CompuCellSetup


sim,simthread = CompuCellSetup.getCoreSimulationObjects()
        
# add extra attributes here
        
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()


from Blastocyst2DSteppables import MitosisSteppable
MitosisSteppableInstance=MitosisSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(MitosisSteppableInstance)
        

from Blastocyst2DSteppables import Original_Config
instanceOfOriginal_Config=Original_Config(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(instanceOfOriginal_Config)


from Blastocyst2DSteppables import Fluid_Growth
instanceOfFluid_Growth=Fluid_Growth(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(instanceOfFluid_Growth)


from Blastocyst2DSteppables import GRN
instanceOfGRN=GRN(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(instanceOfGRN)


from Blastocyst2DSteppables import GRN_Field
instanceOfGRN_Field=GRN_Field(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(instanceOfGRN_Field)


from Blastocyst2DSteppables import Boundary
instanceOfBoundary=Boundary(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(instanceOfBoundary)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        