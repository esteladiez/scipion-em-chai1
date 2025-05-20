# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     you (you@yourinstitution.email)
# *
# * your institution
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
Describe your python module here:
This module will provide the traditional Hello world example
"""
import os
import subprocess
import pwem.objects as emobj
import requests
from pyworkflow.protocol.params import EnumParam, StringParam, FileParam
from pwem.protocols import EMProtocol
from collections import defaultdict
import zipfile
import json
from .. import Plugin


class Chai1Protocol(EMProtocol):
    """Protocol to run Chai-1."""

    RUN_LOCALLY=0
    RUN_SERVER=1

    IMPORT_FASTA=0
    IMPORT_PDBID=1

    _label = 'Chai-1'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)



    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('runserver', EnumParam,
                      default=self.RUN_LOCALLY,
                      label='Run Format',
                      choices=['Run locally installed Chai-1',
                               'Import predictions run on Chai-1 Server'],
                      display=EnumParam.DISPLAY_HLIST)
        form.addParam('source', EnumParam,
                      condition='runserver==%d'%self.RUN_LOCALLY,
                      label='Input file',
                      choices=['Fasta File',
                               'PDB ID'],
                      default=self.IMPORT_PDBID)
        form.addParam('PDBid', StringParam,
                       condition='source==%d' % self.IMPORT_PDBID,
                       label="Protein Data Bank(PDB) NAME/ID: ", allowsNull=True,
                       help='Write a PDB ID (4-character alphanumeric'
                            'characters; examples: 7PZB, 2HBS).\n You can obtain this '
                            'information at https://www.rcsb.org/')
        form.addParam('FASTA', FileParam,
                      condition='source==%d' % self.IMPORT_FASTA,
                      label='Introduce Fasta File')

        form.addParam('serverfile', FileParam,
                      condition='runserver==%d' % self.RUN_SERVER,
                      label='Import File/Folder',
                      help='Import file/folder with the predictions downloaded '
                           'from Chai-1 server')

    def _insertAllSteps(self):
        runserver=self.runserver.get()
        source=self.source.get()
        if runserver == self.RUN_LOCALLY:
            if source==self.IMPORT_FASTA:
                 self._insertFunctionStep('_predictStructure')
            elif source==self.IMPORT_PDBID:
               self._insertFunctionStep('_downloadFastaFile')
               self._insertFunctionStep('_predictStructure')
        elif runserver == self.RUN_SERVER:
            serverfile=self.serverfile.get()
            self._insertFunctionStep('_getModelFromCIF', serverfile)

    def _getModelFromCIF(self, serverfile):
        temp_Dir = self._getExtraPath()  # Ruta temporal donde se extraerán los archivos
        os.makedirs(temp_Dir, exist_ok=True)  # Asegúrate de que el directorio existe
        outFileNames = []
        try:
            # ABRIR ARCHIVO .ZIP Y OBTENER LA INFORMACIÓN QUE NECESITAMOS.
            with zipfile.ZipFile(serverfile, 'r') as zip_ref:
                # EXTRAER  ARCHIVOS .CIF DEL .ZIP FILE

                cif_files = [file for file in zip_ref.namelist() if file.endswith('.cif')]

                if not cif_files:
                    raise Exception(f"No se encontraron archivos .cif en el archivo ZIP: {serverfile}")

                # Extraer únicamente los archivos .cif al directorio temporal
                extracted_files = []
                for cif_file in cif_files:
                    extracted_path = zip_ref.extract(cif_file, temp_Dir)
                    extracted_files.append(extracted_path)
                    print(f"Archivo extraído: {extracted_path}")

                for i, cif_file in enumerate(extracted_files, start=1):
                    amino_acids, average_lddt_values = self.process_cif_file(cif_file)
                    with open(cif_file, 'r') as f:
                        lines = f.readlines()

                    # Modificar las líneas que comienzan con 'ATOM'
                    atom_index = 0
                    new_lines = []
                    for line in lines:
                        if line.startswith('ATOM '):
                            parts = line.split()
                            number = int(parts[8])
                            if number == atom_index + 1:
                                atom_index = atom_index
                            else:
                                atom_index = atom_index + 1
                            parts[
                                13] = f"{average_lddt_values[atom_index]:.10f}"  # Cambiar el valor 1.00 por average_lddt_values
                            new_line = ' '.join(parts) + '\n'
                            new_lines.append(new_line)
                        else:
                            new_lines.append(line)

                    # Escribir el contenido modificado de nuevo en el archivo
                    with open(cif_file, 'w') as f:
                        f.writelines(new_lines)

                    with open(cif_file, 'a') as f:
                        f.write("\n# LDDT values\n")
                        f.write("loop_\n")
                        f.write("_scipion_attributes.name\n")
                        f.write("_scipion_attributes.recipient\n")
                        f.write("scipion_attributes.specifier\n")
                        f.write("scipion_attributes.value\n")
                        cont = 1
                        for index, (amino_acid, lddt_value) in enumerate(zip(amino_acids, average_lddt_values),
                                                                         start=1):
                            f.write(f'LDDT residues: {cont} {amino_acid} {lddt_value}\n')
                            cont = cont + 1

                lddt_files = []
                for i, cif_files in enumerate(extracted_files, start=1):
                    amino_acids, average_lddt_values = self.process_cif_file(cif_files)
                    file_name = self._getPath(f'extra/lddt_{i - 1}.txt')
                    cont = 1
                    with open(file_name, 'w') as f:
                        for amino_acid, lddt_value in zip(amino_acids, average_lddt_values):
                            f.write(f'LDDT  {amino_acid}  {lddt_value:.8f}\n')
                            cont = cont + 1
                        print(f'Archivo {file_name} guardado correctamente.')

                    lddt_files.append(file_name)

                # OBTENER RANKING SCORES
                json_files = [file for file in zip_ref.namelist() if file.endswith('.json') and 'summary' in file]

                if not json_files:
                    raise Exception(f"Not .json files found in the ZIP file: {serverfile}")
                ranking_scores = []
                for json_file in json_files:
                    # Leer el contenido del archivo JSON
                    with zip_ref.open(json_file) as file:
                        data = json.load(file)  # Cargar el JSON como diccionario
                        if 'ranking_score' in data:
                            ranking_scores.append(data['ranking_score'])

                # Crear la información para el summary
                self.clusteringSummary = String()
                clusteringSummary = ''
                for i in range(len(ranking_scores)):
                    msg = '| MODEL ' + str(i) + ' has a ranking score of ' + str(ranking_scores[i])
                    clusteringSummary += msg
                self.clusteringSummary.set(clusteringSummary)
                self._store(self.clusteringSummary)
                print(f"{clusteringSummary}")

            files = glob.glob(os.path.join(temp_Dir, "*.cif"))

            if not files:
                raise Exception("No atomic model selected. No valid .cif files were extracted.")

        except zipfile.BadZipFile:
            raise Exception(f"El archivo proporcionado no es un ZIP válido: {serverfile}")
        except FileNotFoundError:
            raise Exception(f"No se encontró el archivo ZIP: {serverfile}")
        except Exception as e:
            raise Exception(f"Error inesperado: {e}")

    def _downloadFastaFile(self):
        pdb_id = self.PDBid.get()
        print(pdb_id)
        url = f'https://www.rcsb.org/fasta/entry/{pdb_id}'
        fasta_filename = f'{pdb_id}.fasta'
        output_path = os.path.join(self._getTmpPath(), fasta_filename)
        response = requests.get(url)
        if response.status_code == 200:
            with open(output_path, 'w') as f:
                f.write(response.text)
            self._reformatFastaFile(output_path)
            self.info(f"FASTA file downloaded and reformatted to: {output_path}")
            return output_path
        else:
            raise Exception(f"No se pudo descargar el archivo FASTA para {pdb_id}. "
                             f"Código de estado: {response.status_code}")

    def _reformatFastaFile(self, file_path):
        reformatted_lines = []
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('>'):
                    # Extract the protein name
                    parts = line.split('|')
                    if len(parts) > 2:
                        protein_name = parts[2].strip()
                    else:
                        protein_name = "Unknown protein"
                    reformatted_lines.append(f">protein|{protein_name}\n")
                else:
                    reformatted_lines.append(line)

        # Overwrite the file with reformatted content
        with open(file_path, 'w') as f:
            f.writelines(reformatted_lines)

    def _predictStructure(self):
        pdb_id = self.PDBid.get()
        fasta_path = os.path.abspath(os.path.join(self._getTmpPath(), f"{pdb_id}.fasta"))
        output_dir = os.path.abspath(self._getExtraPath())  # Directorio absoluto para guardar el resultado
        os.makedirs(output_dir, exist_ok=True)

        # Verificar si el archivo FASTA existe
        if not os.path.exists(fasta_path):
            raise Exception(f"No se encontró el archivo FASTA en {fasta_path}")

        # Ejecutar el comando de predicción
        predict_args = ["fold", fasta_path, output_dir]
        self.info(f"Ejecutando chai-lab con los siguientes argumentos: {predict_args}")
        Plugin.runChai(self, "chai-lab", predict_args, cwd=output_dir)

        # Verificar los archivos generados en el directorio de salida
        generated_files = os.listdir(output_dir)
        self.info(f"Archivos generados en {output_dir}: {generated_files}")

        # Buscar archivos .cif generados en el directorio
        cif_files = [file for file in generated_files if file.endswith(".cif")]

        if cif_files:
            # Llamar a createOutputStep con los archivos .cif
            atomStructPaths = [os.path.join(output_dir, file) for file in cif_files]
            print(atomStructPaths)
            self.createOutputStep(atomStructPaths)
            self.info(f"Estructuras predichas guardadas: {atomStructPaths}")
        else:
            raise Exception(f"No se encontró el archivo de estructura predicha en {output_dir}")
        for i, cif_file in enumerate(atomStructPaths, start=1):
            amino_acids, average_lddt_values = self.process_cif_file(cif_file)
            with open(cif_file, 'r') as f:
                lines = f.readlines()

            # Modificar las líneas que comienzan con 'ATOM'
            atom_index = 0
            new_lines = []
            for line in lines:
                if line.startswith('ATOM '):
                    parts = line.split()
                    number = int(parts[7])
                    if number == atom_index + 1:
                        atom_index = atom_index
                    else:
                        atom_index = atom_index + 1
                    parts[
                        14] = f"{average_lddt_values[atom_index]:.10f}"  # Cambiar el valor 1.00 por average_lddt_values
                    new_line = ' '.join(parts) + '\n'
                    new_lines.append(new_line)
                else:
                    new_lines.append(line)

            # Escribir el contenido modificado de nuevo en el archivo
            with open(cif_file, 'w') as f:
                f.writelines(new_lines)

            with open(cif_file, 'a') as f:
                f.write("\n# LDDT values\n")
                f.write("loop_\n")
                f.write("_scipion_attributes.name\n")
                f.write("_scipion_attributes.recipient\n")
                f.write("scipion_attributes.specifier\n")
                f.write("scipion_attributes.value\n")
                cont = 1
                for index, (amino_acid, lddt_value) in enumerate(zip(amino_acids, average_lddt_values),
                                                                 start=1):
                    f.write(f'LDDT residues: {cont} {amino_acid} {lddt_value}\n')
                    cont = cont + 1

        lddt_files = []
        for i, cif_files in enumerate(atomStructPaths, start=1):
            amino_acids, average_lddt_values = self.process_cif_file(cif_files)
            file_name = self._getPath(f'extra/lddt_{i - 1}.txt')
            cont = 1
            with open(file_name, 'w') as f:
                for amino_acid, lddt_value in zip(amino_acids, average_lddt_values):
                    f.write(f'LDDT  {amino_acid}  {lddt_value:.8f}\n')
                    cont = cont + 1
                print(f'Archivo {file_name} guardado correctamente.')

            lddt_files.append(file_name)

    def process_cif_file(self, cif_files):
        amino_acid_values = defaultdict(list) # Diccionario para agrupar valores por aminoácido
        amino_acid_names=defaultdict(list)
        cont=1
        with open(cif_files, 'r') as cif_file:
            for line in cif_file:
                if line.startswith("ATOM"):
                    parts = line.split()
                    amino_acid = parts[6]  # Nombre del aminoácido
                    number_aa = int(parts[7])
                    value = float(parts[18])
                    if number_aa == cont:
                        amino_acid_values[cont].append(value)
                        amino_acid_names[cont].append(amino_acid)
                    else:
                        cont=cont+1
                        amino_acid_values[cont].append(value)
                        amino_acid_names[cont].append(amino_acid)
        # Calcular el promedio por aminoácido
        amino_acids = []
        average_values = []
        for number_aa in sorted(amino_acid_values.keys()):  # Ordenar por número de aminoácido
            values = amino_acid_values[number_aa]
            names = amino_acid_names[number_aa]

            avg_value = sum(values) / len(values) if values else 0  # Calcular el promedio de los valores pLDDT
            amino_acids.append(f"{names[0]}")  # Nombre del primer aminoácido
            average_values.append(avg_value)  # Agregar el promedio calculado

        return amino_acids, average_values


    def createOutputStep(self, atomStructPaths, paeFns=[]):
            """ Copy the atomic structure and register the output object.
            :param list_string atomStructPath: list of atom struct files to be
                                                saved
            """
            # Preparar el diccionario de salida
            kwargs = {}

            for atomStructPath in atomStructPaths:
                # Verifica si el archivo existe
                if not os.path.exists(atomStructPath):
                    raise Exception(f"Atomic structure not found at *{atomStructPath}*")

                # Procesa los archivos .pdb o .cif
                if atomStructPath.endswith(".pdb") or atomStructPath.endswith(".cif"):
                    # Crear el objeto de estructura atómica
                    pdb = emobj.AtomStruct()
                    pdb.setFileName(atomStructPath)
                    atomStructPath = os.path.basename(atomStructPath)

                    # Si es un archivo .cif, procesa el nombre correctamente
                    if atomStructPath.endswith(".cif"):
                        keyword = atomStructPath.split(".cif")[0].replace(".", "_")
                    else:
                        keyword = atomStructPath.split(".pdb")[0].replace(".", "_")

                    # Verifica si la palabra clave comienza con un número, y agrega un prefijo si es necesario
                    if keyword[0].isdigit():
                        keyword = "AS_" + keyword
                    kwargs[keyword] = pdb

            # Si se tiene algún archivo .json de PAE, también se registra
            for paeFn in paeFns:
                paeObject = PAE(filename=paeFn)
                paeFn = os.path.basename(paeFn)
                if paeFn.endswith(".json"):
                    keyword = paeFn.split(".json")[0].replace(".", "_")
                else:
                    keyword = paeFn.split(".jsn")[0].replace(".", "_")
                if keyword[0].isdigit():
                    keyword = "PAE_" + keyword
                kwargs[keyword] = paeObject

            # Definir los objetos de salida
            self._defineOutputs(**kwargs)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []


    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        if self.isFinished():
            summary.append(f"This protocol has printed *{self.message}* {self.times} times.")
        return summary

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append(f"{self.message} has been printed in this run {self.times} times.")
            if self.previousCount.hasPointer():
                methods.append("Accumulated count from previous runs were %i."
                               " In total, %s messages has been printed."
                               % (self.previousCount, self.count))
        return methods
