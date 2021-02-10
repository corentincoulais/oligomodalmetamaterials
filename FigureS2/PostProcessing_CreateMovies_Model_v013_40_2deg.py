
import numpy as np
import os
import time

from numpy import sin, cos, tan, sqrt, pi

Model_version = 13
Subversion = 40
#MATID = 1



for MATID in range(18):

    Mdb()

    # -*- coding: mbcs -*-
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


    FrameRate=50

    versioning_name1 = '_v0'+str(Model_version)+'_0'+str(Subversion)
    versioning_name2 = versioning_name1+'-Mat-'+str(MATID)

    file_location =r"F:\15_Viscoelastic_Metamaterial\AbaqusNew\Model_v013\batch_v13_40_2deg_Visco-Agilus_improved_norelax19mm"
    ODBfile_name = file_location+'\\Job-NLViscoComp'+versioning_name2+'.odb'
    os.chdir( file_location)


    output_Movie_name = file_location+'\\MOVIE_Job-NLViscoComp'+versioning_name2

    ########################
    # Save displacement data
    ########################

    o1 = session.openOdb(    name=ODBfile_name)
    odb = session.odbs[ODBfile_name]
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.viewports['Viewport: 1'].setValues(displayedObject=odb)


    session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
        visibleEdges=FEATURE)
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
        CONTOURS_ON_DEF, ))
        
        
    session.animationController.setValues(animationType=TIME_HISTORY, viewports=(
        'Viewport: 1', ))
    session.animationController.play(duration=UNLIMITED)
    session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
    session.imageAnimationOptions.setValues(vpDecorations=OFF, vpBackground=OFF, 
        compass=OFF, timeScale=1, frameRate=FrameRate)

        
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
        scratchCoordSystemDisplay=OFF, pointElements=OFF)
    session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(triad=OFF, 
        legend=OFF, title=OFF, state=OFF, annotations=OFF, compass=OFF)
        
    session.writeImageAnimation(
        fileName=output_Movie_name, 
        format=QUICKTIME, canvasObjects=(session.viewports['Viewport: 1'], ))

    session.odbs[ODBfile_name].close()
    
    time.sleep(3)
