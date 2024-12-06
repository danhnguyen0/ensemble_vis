import sys
import vtk
import math

def main(filename):
    # Colors and Renderer
    colors = vtk.vtkNamedColors()
    
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d("SlateGray"))
    
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(640, 480)
    
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderWindow)

    # Read PDB file
    pdb = vtk.vtkPDBReader()
    pdb.SetFileName(filename)
    pdb.SetHBScale(1.0)
    pdb.SetBScale(1.0)
    pdb.Update()
    print("# of atoms is:", pdb.GetNumberOfAtoms())

    # Calculate resolution based on atom count
    resolution = math.sqrt(300000.0 / pdb.GetNumberOfAtoms())
    resolution = min(max(resolution, 4), 20)
    print("Resolution is:", resolution)

    # Sphere source for atoms
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(0, 0, 0)
    sphere.SetRadius(1)
    sphere.SetThetaResolution(int(resolution))
    sphere.SetPhiResolution(int(resolution))

    # Glyph for atom rendering
    glyph = vtk.vtkGlyph3D()
    glyph.SetInputConnection(pdb.GetOutputPort())
    glyph.SetOrient(1)
    glyph.SetColorMode(1)
    glyph.SetScaleMode(2)  # Scale by scalar value
    glyph.SetScaleFactor(0.25)
    glyph.SetSourceConnection(sphere.GetOutputPort())

    atomMapper = vtk.vtkPolyDataMapper()
    atomMapper.SetInputConnection(glyph.GetOutputPort())
    atomMapper.UseLookupTableScalarRangeOff()
    atomMapper.ScalarVisibilityOn()
    atomMapper.SetScalarModeToDefault()

    atom = vtk.vtkLODActor()
    atom.SetMapper(atomMapper)
    atom.GetProperty().SetRepresentationToSurface()
    atom.GetProperty().SetInterpolationToGouraud()
    atom.GetProperty().SetAmbient(0.1)
    atom.GetProperty().SetDiffuse(0.7)
    atom.GetProperty().SetSpecular(0.5)
    atom.GetProperty().SetSpecularPower(80)
    atom.GetProperty().SetSpecularColor(colors.GetColor3d("White"))
    atom.SetNumberOfCloudPoints(30000)

    renderer.AddActor(atom)

    # Tube filter for bonds
    tube = vtk.vtkTubeFilter()
    tube.SetInputConnection(pdb.GetOutputPort())
    tube.SetNumberOfSides(int(resolution))
    tube.CappingOff()
    tube.SetRadius(0.2)
    tube.SetVaryRadius(0)
    tube.SetRadiusFactor(10)

    bondMapper = vtk.vtkPolyDataMapper()
    bondMapper.SetInputConnection(tube.GetOutputPort())
    bondMapper.UseLookupTableScalarRangeOff()
    bondMapper.ScalarVisibilityOff()
    bondMapper.SetScalarModeToDefault()

    bond = vtk.vtkLODActor()
    bond.SetMapper(bondMapper)
    bond.GetProperty().SetRepresentationToSurface()
    bond.GetProperty().SetInterpolationToGouraud()
    bond.GetProperty().SetAmbient(0.1)
    bond.GetProperty().SetDiffuse(0.7)
    bond.GetProperty().SetSpecular(0.5)
    bond.GetProperty().SetSpecularPower(80)
    bond.GetProperty().SetSpecularColor(colors.GetColor3d("White"))

    renderer.AddActor(bond)

    # Render and start interaction
    renderWindow.SetWindowName("ReadPDB")
    renderWindow.Render()
    interactor.Initialize()
    interactor.Start()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python", sys.argv[0], "Filename(.pdb) e.g. PED00020e001.pdb")
        sys.exit(1)
    main(sys.argv[1])
