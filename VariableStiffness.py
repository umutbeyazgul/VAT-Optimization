#
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from math import *
import numpy as np
import time
from scipy.optimize import differential_evolution
initialtime = time.time()
# plate sizes x, y
xlength=500
ylength=500
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-xlength*0.5, -ylength*0.5), 
    point2=(xlength*0.5, ylength*0.5))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseShell(sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

# meshing
boyut=25 #element size
mdb.models['Model-1'].parts['Part-1'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=boyut)
mdb.models['Model-1'].parts['Part-1'].generateMesh()

#Material mechanical properties Mpa
mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Elastic(table=((127700.0, 7400.0, 
    0.33, 6900.0, 6900.0, 4300.0), ), type=LAMINA)
t=0.2 #thickness mm
numberofplies=16
#creating sets for elements
ElemanSayisi=len(mdb.models['Model-1'].parts['Part-1'].elements)
for i in range(ElemanSayisi):
    mdb.models['Model-1'].parts['Part-1'].Set(elements=
        mdb.models['Model-1'].parts['Part-1'].elements[i:i+1], name = "Set-%s" % (i+1))
#
# Ply Properties. CompositeLayup
#CompositeLayup
for i in range(ElemanSayisi):
    mdb.models['Model-1'].parts['Part-1'].CompositeLayup(description='', 
        elementType=SHELL, name='CompositeLayup-%s'%(i+1), offsetType=MIDDLE_SURFACE, 
        symmetric=False, thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1'].compositeLayups['CompositeLayup-%s'%(i+1)].Section(
        integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, 
        temperature=GRADIENT, thicknessType=UNIFORM, useDensity=OFF)
    mdb.models['Model-1'].parts['Part-1'].compositeLayups['CompositeLayup-%s'%(i+1)].ReferenceOrientation(
        additionalRotationType=ROTATION_NONE, angle=0.0, axis=AXIS_3, fieldName='', 
        localCsys=None, orientationType=GLOBAL)
###########################################################################################
def Buckling(T):
    for i in range(ElemanSayisi):
        mdb.models['Model-1'].parts['Part-1'].compositeLayups['CompositeLayup-%s'%(i+1)].deletePlies() # delete previous plies
    #Layer Sequence ve orientation angles:
    #Single angles:[a] constant orientation angle; double angles:[a, b] variable orientation angle
    plylist1sthalf=[[T[0], T[1]],[-T[0], -T[1]],[T[2], T[3]],[-T[2], -T[3]],[T[4], T[5]],[-T[4], -T[5]],[T[6], T[7]],[-T[6], -T[7]]]
    PlylistSymm=plylist1sthalf+plylist1sthalf[::-1] #simetrik dizilim
    #Plies and orientation distribution
    plynum=range(len(PlylistSymm))
    for phi, j in zip(PlylistSymm, plynum):
        for i in range(ElemanSayisi):
            if len(phi)==2:
                coord0=mdb.models['Model-1'].parts['Part-1'].sets["Set-%s" % (i+1)].nodes[0].coordinates #elementwise node numbering 0,1,2,3
                coord1=mdb.models['Model-1'].parts['Part-1'].sets["Set-%s" % (i+1)].nodes[1].coordinates
                coord2=mdb.models['Model-1'].parts['Part-1'].sets["Set-%s" % (i+1)].nodes[2].coordinates
                coord3=mdb.models['Model-1'].parts['Part-1'].sets["Set-%s" % (i+1)].nodes[3].coordinates
                coordlist=[coord0,coord1,coord2,coord3]
                elemcoord=[sum(x)/len(x) for x in zip(*coordlist)] #element center
                theta=2*abs(elemcoord[0])*(phi[1]-phi[0])/xlength+phi[0] #linear angle orientation
            else:
                theta=phi[0]
            #
            mdb.models['Model-1'].parts['Part-1'].compositeLayups['CompositeLayup-%s'%(i+1)].CompositePly(
                additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
                , axis=AXIS_3, material='Material-1', numIntPoints=3, orientationType=
                SPECIFY_ORIENT, orientationValue=theta, plyName='Ply-%s'%(j+1), region=
            mdb.models['Model-1'].parts['Part-1'].sets["Set-%s" % (i+1)], suppressed=False, 
                thickness=t, thicknessType=SPECIFY_THICKNESS)
    #
    #Assembly
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
        part=mdb.models['Model-1'].parts['Part-1'])
    #############
    mdb.models['Model-1'].BuckleStep(maxIterations=300, name='Step-1', numEigen=1, 
        previous='Initial', vectors=2)
    mdb.models['Model-1'].rootAssembly.Surface(name='Surf-1', side1Edges=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.getSequenceFromMask(
        ('[#4 ]', ), ))
    mdb.models['Model-1'].ShellEdgeLoad(createStepName='Step-1', distributionType=
        UNIFORM, field='', localCsys=None, magnitude=1.0, name='Load-1', region=
        mdb.models['Model-1'].rootAssembly.surfaces['Surf-1'])
    mdb.models['Model-1'].rootAssembly.Set(edges=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.getSequenceFromMask(
        ('[#e ]', ), ), name='Set-1')
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
        region=mdb.models['Model-1'].rootAssembly.sets['Set-1'], u1=UNSET, u2=UNSET
        , u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mdb.models['Model-1'].rootAssembly.Set(edges=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.getSequenceFromMask(
        ('[#1 ]', ), ), name='Set-2')
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
        region=mdb.models['Model-1'].rootAssembly.sets['Set-2'], u1=SET, u2=SET, 
        u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name='Job-1', nodalOutputPrecision=SINGLE, 
        numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
        ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
    # Sonuc
    mdb.jobs['Job-1'].waitForCompletion()
    odbAsVariable=session.openOdb(name='Job-1.odb') 
    Frames=odbAsVariable.steps['Step-1'].frames 
    FrameInfo=Frames[1].description 
    infolist=FrameInfo.split() 
    BucklingLoad=-1*float(infolist[-1]) 
    f.write(str(plylist1sthalf)+"\t"+str(BucklingLoad)+"\n")
    return BucklingLoad
#Optimization
f= open("BurkulmaYukleri.txt","w")
bounds=[] #lower and upper bounds
for i in range(int(numberofplies/2)): #number of variable
    bounds = bounds+[(-90,90)]
result = differential_evolution(Buckling, bounds, popsize=5, tol=0.01) #
print(str(result.x), str(result.fun),"\n")
f.write("Result:"+str(result.x)+"\t"+str(result.fun)+"\n")
elapsed = time.time() - initialtime
f.write("Elapsed Time:"+"\t"+str(elapsed)+"\n")
f.close()
print(elapsed)