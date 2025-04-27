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
import requests
from pyworkflow.protocol.params import EnumParam, StringParam, FileParam
from pwem.protocols import EMProtocol
from .. import Plugin


class Chai1Protocol(EMProtocol):
    """Protocol to run Chai-1."""

    IMPORT_FASTA=0
    IMPORT_PDBID=1

    _label = 'Chai-1'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)



    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('source', EnumParam,
                      choices=['Fasta File',
                               'PDB ID'],
                      display=EnumParam.DISPLAY_HLIST,)
        form.addParam('PDBid', StringParam,
                      condition='source==%d' % self.IMPORT_PDBID,
                       label="Protein Data Bank(PDB) NAME/ID: ", allowsNull=True,
                       help='Write a PDB ID (4-character alphanumeric'
                            'characters; examples: 7PZB, 2HBS).\n You can obtain this '
                            'information at https://www.rcsb.org/')
        form.addParam('FASTA', FileParam,
                      condition='source==%d' % self.IMPORT_FASTA,
                      label='Introduce Fasta File')

    def _insertAllSteps(self):
        source=self.source.get()
        if source==self.IMPORT_FASTA:
            self._insertFunctionStep('_predictStructure')
        elif source==self.IMPORT_PDBID:
           self._insertFunctionStep('_downloadFastaFile')
           self._insertFunctionStep('_predictStructure')

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
        fasta_path = os.path.join(self._getTmpPath(), f"{pdb_id}.fasta")
        output_dir = self._getExtraPath()  # Directorio para guardar resultado
        os.makedirs(output_dir, exist_ok=True)

        # Ejecutar el comando de predicción
        predict_args = ["fold", fasta_path, output_dir]
        Plugin.runChai(self, "chai-lab", predict_args, cwd=output_dir)

        # Asumimos que genera un archivo .pdb con nombre <PDB_ID>_predicted.pdb
        predicted_file = os.path.join(output_dir, f"{pdb_id}_predicted.pdb")
        if os.path.exists(predicted_file):
            self.info(f"Estructura predicha guardada en: {predicted_file}")
        else:
            raise Exception(f"No se encontró el archivo de estructura predicha en {output_dir}")
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
