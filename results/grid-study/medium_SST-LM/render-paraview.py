"""Render the reconstructed OpenFOAM velocity field with ParaView 6."""

from pathlib import Path
import sys

from paraview import simple


case = Path(sys.argv[1]).resolve()
iteration = int(sys.argv[2])
output = Path(sys.argv[3]).resolve()
output.mkdir(parents=True, exist_ok=True)
marker = case / 'baseline.foam'
marker.touch(exist_ok=True)

reader = simple.OpenFOAMReader(registrationName='Baseline solution', FileName=str(marker))
reader.MeshRegions = ['internalMesh']
reader.CellArrays = ['U']
reader.UpdatePipeline(iteration)
if iteration not in reader.TimestepValues:
    raise ValueError(f'Iteration {iteration} is missing from {marker}.')
bounds = reader.GetDataInformation().GetBounds()
midspan = 0.5 * (bounds[4] + bounds[5])

plane = simple.Slice(registrationName='Midspan', Input=reader)
plane.SliceType = 'Plane'
plane.SliceType.Origin = [0, 0, midspan]
plane.SliceType.Normal = [0, 0, 1]
plane.UpdatePipeline(iteration)
if not plane.GetDataInformation().GetNumberOfCells():
    raise RuntimeError('ParaView found no cells on the midspan plane.')

view = simple.GetActiveViewOrCreate('RenderView')
view.InteractionMode = '2D'
view.CameraParallelProjection = 1
view.ViewTime = iteration
view.UseColorPaletteForBackground = 0
view.Background = [0.96, 0.97, 0.98]
view.OrientationAxesVisibility = 0
simple.GetAnimationScene().AnimationTime = iteration
simple.GetLayout().SetSize(1400, 760)

display = simple.Show(plane, view)
display.Representation = 'Surface'
simple.ColorBy(display, ('CELLS', 'U', 'Magnitude'))
colors = simple.GetColorTransferFunction('U')
colors.ApplyPreset('Viridis (matplotlib)', True)
colors.RescaleTransferFunction(0.0, 22.0)
colors.AutomaticRescaleRangeMode = 'Never'
display.SetScalarBarVisibility(view, True)
legend = simple.GetScalarBar(colors, view)
legend.Title = 'Velocity magnitude'
legend.ComponentTitle = 'm/s'
legend.Orientation = 'Horizontal'
legend.WindowLocation = 'Lower Center'
legend.ScalarBarLength = 0.39
legend.TitleColor = [1.0, 1.0, 1.0]
legend.LabelColor = legend.TitleColor
legend.TitleFontSize = 17
legend.LabelFontSize = 15

label = simple.Text()
label_display = simple.Show(label, view)
label_display.WindowLocation = 'Upper Left Corner'
label_display.Color = [1.0, 1.0, 1.0]
label_display.FontSize = 22

for name, description, x, y, scale in (
    ('paraview-velocity.png', 'Airfoil and wake', 0.120, 0.008, 0.104),
    ('paraview-leading-edge.png', 'Leading edge', 0.022, 0.006, 0.022),
):
    label.Text = f'NACA 0012  |  SST-LM  |  Iteration {iteration}\n{description} at midspan'
    simple.Render(view)
    view.CameraFocalPoint = [x, y, midspan]
    view.CameraPosition = [x, y, midspan + 1.0]
    view.CameraViewUp = [0, 1, 0]
    view.CameraParallelScale = scale
    simple.SaveScreenshot(str(output / name), view, ImageResolution=[1400, 760])

simple.SaveState(str(output / 'paraview-state.pvsm'))
print(f'ParaView rendered iteration {iteration} at z = {midspan:g} m')
