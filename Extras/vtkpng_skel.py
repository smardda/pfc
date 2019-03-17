#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
jelotvtk = LegacyVTKReader(FileNames=['/home/warter/jetaux/Data/VTK/jelot.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1004, 697]

# get color transfer function/color map for 'Body'
bodyLUT = GetColorTransferFunction('Body')

# show data in view
jelotvtkDisplay = Show(jelotvtk, renderView1)
# trace defaults for the display properties.
jelotvtkDisplay.ColorArrayName = ['CELLS', 'Body']
jelotvtkDisplay.LookupTable = bodyLUT
jelotvtkDisplay.ScalarOpacityUnitDistance = 92.99815822551636

# reset view to fit data
renderView1.ResetCamera()

## show color bar/color legend
#jelotvtkDisplay.SetScalarBarVisibility(renderView1, False)

# get opacity transfer function/opacity map for 'Body'
bodyPWF = GetOpacityTransferFunction('Body')

# turn off scalar coloring
ColorBy(jelotvtkDisplay, None)

# Properties modified on jelotvtkDisplay
jelotvtkDisplay.Opacity = 0.3

# set active source
SetActiveSource(jelotvtk)

# create a new 'Legacy VTK Reader'
partvtk = LegacyVTKReader(FileNames=['PART.vtk'])

# show data in view
partvtkDisplay = Show(partvtk, renderView1)
# trace defaults for the display properties.
partvtkDisplay.ColorArrayName = ['CELLS', 'Body']
partvtkDisplay.LookupTable = bodyLUT
partvtkDisplay.ScalarOpacityUnitDistance = 66.20098695520572

# rescale color and/or opacity maps used to exactly fit the current data range
partvtkDisplay.RescaleTransferFunctionToDataRange(False)

# show color bar/color legend
partvtkDisplay.SetScalarBarVisibility(renderView1, True)

# current camera placement for renderView1
renderView1.CameraPosition = [629.125355665469, -10299.351248630117, 5047.820984626171]
renderView1.CameraFocalPoint = [-234.1169810554541, 4094.2957792777675, -2685.9298473469676]
renderView1.CameraViewUp = [-0.04647285132183669, 0.47063293457535416, 0.8811043723549492]
renderView1.CameraParallelScale = 6200.376426408011

# save screenshot
SaveScreenshot('PART.png', magnification=1, quality=100, view=renderView1)


#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

quit()
