import os

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol import params
from ..protocols.protocol_chai1 import Chai1Protocol
from pwem.viewers import ChimeraAttributeViewer
import matplotlib.pyplot as plt


from pwem.viewers.viewer_chimera import (Chimera,
                                         sessionFile)
from chimera.objects import PAE

class ChimeraChaiViewer(ChimeraAttributeViewer):
    _label = 'viewer chai1'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [Chai1Protocol]
    _viewerOptions = ['ChimeraX']
    _viewmodeloptions = ['MODEL 0', 'MODEL 1', 'MODEL 2', 'MODEL 3', 'MODEL 4']

    def _defineParams(self, form):
        form.addSection(label='Visualization of model in ChimeraX')
        group = form.addGroup('AtomStruct General Viewer')
        group.addParam('displaySoftware', params.EnumParam,
                       choices=self._viewerOptions, default=0,
                       label='Display AtomStruct with: ',
                       help='Display the AtomStruct object with which software.\nAvailable: PyMol, ChimeraX')
        form.addSection(label='Visualization of LDDT information')
        group = form.addGroup('Display LDDT information')
        group.addParam('display', params.BooleanParam, default=False,
                       labe='Do you want to display LDDT information?')
        group.addParam('model', params.EnumParam,
                       label='Which model information do you want to diplay?',
                       choices=self._viewmodeloptions, default=0, condition='display')
        form.addParam('viewConservation', params.LabelParam,
                      label='Display conservation over sequence: ',
                      help='Display a graph with the values of the selected attribute over the sequence.',
                      condition='display')
        form.addParam('information', params.IntParam, default=1,
                      label='Aminoacid Podition:',
                      help='Obtain information about the aminoacid that corresponds to each position',
                      condition='display')
        form.addParam('name', params.LabelParam,
                      label='The name of the aminoacid located at the position specified is: ', condition='display')


def _getVisualizeDict(self):
    return {
        'displaySoftware': self._viewAtomStruct,
        'viewConservation': self._showlddt,
        'name': self._showaminoacid,
    }


def _viewAtomStruct(self, e=None):
    if self.displaySoftware.get() == 0:
        return self._visualize(self.getAtomStruct())


def getAtomStruct(self):
    obj = self.protocol
    return obj


def _showlddt(self, paramName=None):
    obj = self.protocol._outputs
    fileNames = []

    # Obtener las rutas absolutas de los archivos generados
    for output in obj:
        file_path = os.path.abspath(eval(f'self.protocol.{output}.getFileName()'))
        fileNames.append(file_path)

    if not fileNames:
        print("No se encontraron archivos de salida.")
        return

    # Seleccionar el archivo del modelo especificado
    model_index = self.model.get()
    if model_index >= len(fileNames):
        print(f"El modelo seleccionado ({model_index}) no existe.")
        return

    file_path = fileNames[model_index]
    aa = []
    values = []
    number = []
    # Leer el archivo y extraer los datos de LDDT
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('LDDT residues:'):
                parts = line.split()
                number.append(int(parts[2]))
                aa.append(parts[3])  # Aminoácido
                values.append(float(parts[4]))  # Valor LDDT

    if not aa or not values:
        print(f"No se encontraron datos LDDT en {file_path}")
        return

    fig, ax = plt.subplots(figsize=(12, 5))

    # Crear el gráfico de barras
    ax.bar(number, values, color='r', alpha=0.7)

    # Configuración del diseño
    ax.set_xlabel('Aminoacid number', fontsize=12, fontweight='bold')
    ax.set_ylabel('LDDT value', fontsize=12, fontweight='bold')
    ax.set_title(f'LDDT per aminoacid - MODEL: {model_index}', fontsize=14, fontweight='bold')

    # Establecer marcas en el eje X cada 100
    ax.set_xticks(range(0, max(number) + 1, 100))

    # Activar la cuadrícula
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # Ajustar el diseño para evitar solapamiento
    plt.tight_layout()

    # Mostrar la gráfica
    plt.show()


def _showaminoacid(self, event=None):
    obj = self.protocol._outputs
    fileNames = []

    # Obtener las rutas absolutas de los archivos generados
    for output in obj:
        file_path = os.path.abspath(eval(f'self.protocol.{output}.getFileName()'))
        fileNames.append(file_path)

    if not fileNames:
        print("No se encontraron archivos de salida.")
        return

    # Seleccionar el archivo del modelo especificado
    model_index = self.model.get()
    if model_index >= len(fileNames):
        print(f"El modelo seleccionado ({model_index}) no existe.")
        return

    file_path = fileNames[model_index]
    aa = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('LDDT residues:'):
                parts = line.split()
                aa.append(parts[3])  # Aminoácido
    information = self.information.get()
    position = int(information)
    aminoacid = aa[position]
    print(aminoacid)


def _visualize(self, obj, **args):
    # create axis file
    models = 1
    dim = 150
    sampling = 1.
    extraFileName = os.path.abspath(self.protocol._getExtraPath("axis_input.bild"))
    Chimera.createCoordinateAxisFile(dim,
                                     bildFileName=extraFileName,
                                     sampling=sampling)

    fnCmd = self.protocol._getExtraPath("chimera_alphafold.cxc")
    f = open(fnCmd, 'w')
    f.write("open %s\n" % extraFileName)
    models += 1
    f.write("cofr 0,0,0\n")  # set center of coordinates
    # change to workingDir
    # If we do not use cd and the project name has an space
    # the protocol fails even if we pass absolute paths
    f.write('cd %s\n' % os.getcwd())

    # get path to atomstructs
    for output in self.protocol._outputs:
        # if the file is an atomic struct show it in chimera
        fileName = os.path.abspath(eval(f'self.protocol.{output}.getFileName()'))
        if fileName.endswith(".cif") or fileName.endswith(".pdb"):
            f.write("open %s\n" % fileName)
            models += 1
    # if exists upload other results files
    # model_?_unrelaxed.pdb
    # pattern = self.protocol._getExtraPath("results/model_?_unrelaxed.pdb")
    # from glob import glob
    # for model in glob(pattern):
    #    f.write("open %s\n" % model)
    #    f.write(f"hide #{models} models\n")
    #    models +=1
    # set alphafold colormap
    f.write("color bfactor palette alphafold\n")
    f.write("key red:low orange: yellow: cornflowerblue: blue:high\n")
    f.close()
    Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")
