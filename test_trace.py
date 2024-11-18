import sys
import vtk
import math
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QSlider, QLabel, QLineEdit, QPushButton
from PyQt5.QtCore import Qt

class PDBVisualizer(QWidget):
    def __init__(self, filename):
        super().__init__()
        self.filename = filename
        self.max_atoms = 1000  # The maximum number of atoms you want to visualize
        self.slider_value = self.max_atoms  # Initialize slider_value
        self.initUI()  # Initialize the UI components

    def initUI(self):
        self.setWindowTitle('PDB Atom Visualizer')

        # Create label to show current input value
        self.label = QLabel(f'Atoms to visualize: {self.slider_value}', self)

        # Create QLineEdit for user input
        self.input_field = QLineEdit(self)
        self.input_field.setValidator(Qt.QIntValidator(1, self.max_atoms))  # Only allow integers in range
        self.input_field.setText(str(self.slider_value))  # Set initial value
        self.input_field.textChanged.connect(self.updateInput)

        # Create button to trigger atom rendering
        self.render_button = QPushButton('Render Atoms', self)
        self.render_button.clicked.connect(self.render_atoms)

        # Set up layout
        layout = QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.input_field)
        layout.addWidget(self.render_button)
        self.setLayout(layout)

        self.show()
        self.render_atoms()

    def updateSlider(self, value):
        self.slider_value = value
        self.label.setText(f'Atoms to visualize: {value}')
        self.render_atoms()

    def render_atoms(self):
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
        pdb.SetFileName(self.filename)
        pdb.SetHBScale(1.0)
        pdb.SetBScale(1.0)
        pdb.Update()
        print("# of atoms is:", pdb.GetNumberOfAtoms())

        # Get points of the loaded data
        points = pdb.GetOutput().GetPoints()
        total_atoms = points.GetNumberOfPoints()

        # Limit the number of atoms to visualize
        num_atoms_to_render = min(total_atoms, self.slider_value)
        print(f"Rendering {num_atoms_to_render} atoms.")

        # Create a new vtkPoints object to hold the subset of points
        selected_points = vtk.vtkPoints()

        # Add a subset of atoms to the selected_points object
        for i in range(num_atoms_to_render):
            selected_points.InsertNextPoint(points.GetPoint(i))

        # Create a new polydata object with the selected points
        selected_polydata = vtk.vtkPolyData()
        selected_polydata.SetPoints(selected_points)

        # Sphere source for atoms
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(0, 0, 0)
        sphere.SetRadius(1)
        sphere.SetThetaResolution(10)
        sphere.SetPhiResolution(10)

        # Glyph for atom rendering
        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(selected_polydata)
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

        renderer.AddActor(atom)
        
        tube = vtk.vtkTubeFilter()
        tube.SetInputConnection(sphere.GetOutputPort())
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

    app = QApplication(sys.argv)
    visualizer = PDBVisualizer(sys.argv[1])
    sys.exit(app.exec_())
