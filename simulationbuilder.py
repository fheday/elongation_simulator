'''
Gui module to build and save a simulation set
'''

import os
from PyQt5.QtWidgets import QApplication, QWidget, QGridLayout, QLabel, QPushButton,\
    QFileDialog, QSpinBox, QDoubleSpinBox, QCheckBox, QListWidget, QTableWidget,\
    QHeaderView, QTableWidgetItem
from Bio import SeqIO

class SimulationBuilder():
    """
    class to generate GUI.
    """
    def __init__(self):
        self.pre_populate = False
        self.concentrations_file_path = None
        self.genes_dict = {}

        self.app = QApplication([])
        self.window = QWidget()
        grid = QGridLayout()
        self.window.setLayout(grid)


        # self.window.setGeometry(50, 50, 320, 200)
        self.window.setWindowTitle("PyQt5 Example")



        grid.addWidget(QLabel("concentration file: "), 0, 0)
        selected_concentration_label = QLabel("      ")
        selected_concentration_label.setObjectName("selected_concentration_label")
        grid.addWidget(selected_concentration_label, 0, 1)

        concentration_file_open_button = QPushButton('Open file')
        concentration_file_open_button.setToolTip('Open concentration file')
        concentration_file_open_button.clicked.connect(self.open_concentrations_file)
        grid.addWidget(concentration_file_open_button, 0, 2)

        grid.addWidget(QLabel("Initiation rate: "), 1, 0)
        init_rate_spinbox = QDoubleSpinBox()
        init_rate_spinbox.setObjectName("init_rate_spinbox")
        init_rate_spinbox.setRange(0.01, 100.00)
        init_rate_spinbox.setSingleStep(0.01)
        grid.addWidget(init_rate_spinbox, 1, 1)

        grid.addWidget(QLabel("Termination rate: "), 2, 0)
        term_rate_spinbox = QDoubleSpinBox()
        term_rate_spinbox.setObjectName("term_rate_spinbox")
        term_rate_spinbox.setRange(0.01, 100.00)
        term_rate_spinbox.setSingleStep(0.01)
        grid.addWidget(term_rate_spinbox, 2, 1)

        grid.addWidget(QLabel("Pre-populate: "), 3, 0)
        pre_populate_check_box = QCheckBox()
        pre_populate_check_box.setObjectName("pre_populate_check_box")
        grid.addWidget(pre_populate_check_box, 3, 1)

        grid.addWidget(QLabel("Add FASTA file :"), 4, 0)
        add_fasta_file_button = QPushButton("Add file")
        add_fasta_file_button.setToolTip(\
            "Reads FASTA file and make its genes available for simulation")
        add_fasta_file_button.clicked.connect(self.open_fasta_file)
        grid.addWidget(add_fasta_file_button, 4, 1)

        grid.addWidget(QLabel("Loaded FASTA files:"), 5, 0)
        fasta_files_listbox = QListWidget()
        fasta_files_listbox.setObjectName("fasta_files_listbox")
        fasta_files_listbox.clicked.connect(self.onselect_fasta_file)
        grid.addWidget(fasta_files_listbox, 5, 1)

        grid.addWidget(QLabel("Genes:"), 5, 2)
        genes_listbox = QListWidget()
        genes_listbox.setObjectName("genes_listbox")
        grid.addWidget(genes_listbox, 5, 3)

        grid.addWidget(QLabel("Transcript copy number: "), 5, 4)
        gene_copy_number_spinbox = QSpinBox()
        genes_listbox.setSelectionMode(2)
        gene_copy_number_spinbox.setObjectName("gene_copy_number_spinbox")
        gene_copy_number_spinbox.setRange(1, 1000)
        gene_copy_number_spinbox.setSingleStep(1)
        grid.addWidget(gene_copy_number_spinbox, 5, 5)

        add_gene_button = QPushButton("Add gene")
        add_gene_button.clicked.connect(self.add_simulation_entry)
        grid.addWidget(add_gene_button, 6, 1, 1, 4)

        grid.addWidget(QLabel("Added simulations: "), 7, 0)
        added_simulations_listbox = QTableWidget()
        added_simulations_listbox.setObjectName("added_simulations_listbox")
        added_simulations_listbox.setColumnCount(3)
        added_simulations_listbox.setShowGrid(False)
        added_simulations_listbox.setHorizontalHeaderLabels(['Fasta file', 'Gene', 'Copy Number'])
        added_simulations_listbox.resizeColumnsToContents()
        added_simulations_listbox.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        added_simulations_listbox.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        grid.addWidget(added_simulations_listbox, 8, 0, 6, 6)

        remove_entry_button = QPushButton("Remove entry")
        remove_entry_button.clicked.connect(self.remove_simulation_entry)
        grid.addWidget(remove_entry_button, 15, 5)
        self.show()

    def show(self):
        """
        bring GUI window to front and show it.
        """
        self.window.show()
        self.app.exec_()

    def open_concentrations_file(self):
        """
        open dialog to select concentration file.
        """
        concentrations_file, _ = QFileDialog.getOpenFileName(self.window, "Select file",
                                                             os.getcwd(), "CSV File (*.csv )")
        if concentrations_file == '':
            return
        selected_concentration_label = self.window.findChild(QLabel, "selected_concentration_label")
        selected_concentration_label.setText(concentrations_file)


    def open_fasta_file(self):
        """
        Open dialog to select FASTA file with genes to simulate.
        Also populate gene box.
        """
        fasta_file_chooser, _ = QFileDialog\
            .getOpenFileName(self.window, "Select file", os.getcwd(),
                             "txt File (*.txt );;FASTA file (*.fasta);;FNA file (*.fna)")
        if fasta_file_chooser == '' or fasta_file_chooser in self.genes_dict.keys():
            return
        self.genes_dict[fasta_file_chooser] = \
            [record.id for record in SeqIO.parse(fasta_file_chooser, "fasta")]
        #update fasta files listbox
        fasta_files_listbox = self.window.findChild(QListWidget, 'fasta_files_listbox')
        if fasta_files_listbox.count() > 0:
            fasta_files_listbox.clear()
        fasta_files_listbox.addItems(self.genes_dict.keys())
        fasta_files_listbox.setCurrentRow(fasta_files_listbox.count() - 1)
        #update genes listbox
        genes_listbox = self.window.findChild(QListWidget, 'genes_listbox')
        if genes_listbox.count() > 0:
            genes_listbox.clear()
        genes_listbox.addItems(self.genes_dict[fasta_file_chooser])
        genes_listbox.setCurrentRow(genes_listbox.count())

    def onselect_fasta_file(self):
        """
        Update gene box with the genes of the selected fasta file.
        """
        fasta_files_listbox = self.window.findChild(QListWidget, 'fasta_files_listbox')
        value = fasta_files_listbox.currentItem().text()
        #update genes listbox
        genes_listbox = self.window.findChild(QListWidget, 'genes_listbox')
        if genes_listbox.count() > 0:
            genes_listbox.clear()
        genes_listbox.addItems(self.genes_dict[value])
        genes_listbox.setCurrentRow(genes_listbox.count())

    def add_simulation_entry(self):
        """
        Use the selected fasta file, gene and gene copy number to add a new simulation entry.
        If fasta file AND gene is already added, update that entry with new copy number.
        """
        # we need: a gene and copy number.
        genes_listbox = self.window.findChild(QListWidget, 'genes_listbox')
        fasta_files_listbox = self.window.findChild(QListWidget, 'fasta_files_listbox')
        added_simulations_listbox = self.window.findChild(QTableWidget, 'added_simulations_listbox')
        gene_copy_number_spinbox = self.window.findChild(QSpinBox, 'gene_copy_number_spinbox')
        if len(genes_listbox.selectedItems()) == 0 or len(fasta_files_listbox.selectedItems()) != 1:
            return
        selected_fasta_file = fasta_files_listbox.currentItem().text()
        selected_genes = [item.text() for item in genes_listbox.selectedItems()]
        fasta_files_list = [str(added_simulations_listbox.item(i, 0).text())
                            for i in range(added_simulations_listbox.rowCount())]
        genes_list = [str(added_simulations_listbox.item(i, 1).text())
                      for i in range(added_simulations_listbox.rowCount())]

        fasta_files_list_indexes = [i for i, x in enumerate(fasta_files_list)
                                    if x == selected_fasta_file]
        genes_list_indexes = [i for i, x in enumerate(genes_list) if x in selected_genes]
        index = set(fasta_files_list_indexes).intersection(genes_list_indexes)

        if len(index) > 0:
            # update gene copy number
            for i in index:
                added_simulations_listbox.setItem(i,
                                                2,
                                                QTableWidgetItem(
                                                    str(gene_copy_number_spinbox.value())))
        else:
            # new entry
            for selected_gene in selected_genes:
                index = added_simulations_listbox.rowCount()
                added_simulations_listbox.insertRow(index)
                added_simulations_listbox.setItem(index, 0, QTableWidgetItem(selected_fasta_file))
                added_simulations_listbox.setItem(index, 1, QTableWidgetItem(selected_gene))
                added_simulations_listbox.setItem(index, 2,
                                                  QTableWidgetItem(
                                                      str(gene_copy_number_spinbox.value())))
        return

    def remove_simulation_entry(self):
        """
        Remove selected simulation entry.
        """
        added_simulations_listbox = self.window.findChild(QTableWidget, 'added_simulations_listbox')
        index = added_simulations_listbox.currentRow()
        added_simulations_listbox.removeRow(index)


if __name__ == "__main__":
    SimulationBuilder()
