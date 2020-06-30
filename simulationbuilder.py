'''
Gui module to build and save a simulation set
'''

import os
from PyQt5.QtWidgets import QApplication, QWidget, QGridLayout, QLabel, QPushButton,\
    QFileDialog, QSpinBox, QDoubleSpinBox, QCheckBox, QListWidget, QTableWidget,\
    QHeaderView, QTableWidgetItem, QGroupBox, QRadioButton
from Bio import SeqIO

class Gui():
    """
    class to generate GUI.
    """
    def __init__(self):
        self.pre_populate = False
        self.concentrations_file_path = None
        self.genes_dict = {}

        self.app = QApplication([])
        self.window = QWidget()
        self.window.resize(100, 800) #initial window size
        grid = QGridLayout()
        self.window.setLayout(grid)

        self.window.setWindowTitle("Simulation Builder")

        grid.addWidget(QLabel("concentration file: "), 0, 0)
        selected_concentration_label = QLabel("      ")
        selected_concentration_label.setObjectName("selected_concentration_label")
        grid.addWidget(selected_concentration_label, 0, 1)

        concentration_file_open_button = QPushButton('Open file')
        concentration_file_open_button.setToolTip('Open concentration file')
        concentration_file_open_button.clicked.connect(self.open_concentrations_file)
        grid.addWidget(concentration_file_open_button, 0, 2)


        grid.addWidget(QLabel("Pre-populate: "), 1, 0)
        pre_populate_check_box = QCheckBox()
        pre_populate_check_box.setObjectName("pre_populate_check_box")
        grid.addWidget(pre_populate_check_box, 1, 1)

        grid.addWidget(QLabel("Add FASTA file :"), 2, 0)
        add_fasta_file_button = QPushButton("Add file")
        add_fasta_file_button.setToolTip(\
            "Reads FASTA file and make its genes available for simulation")
        add_fasta_file_button.clicked.connect(self.open_fasta_file)
        grid.addWidget(add_fasta_file_button, 2, 1)

        grid.addWidget(QLabel("Loaded FASTA files:"), 3, 0)
        fasta_files_listbox = QListWidget()
        fasta_files_listbox.setObjectName("fasta_files_listbox")
        fasta_files_listbox.clicked.connect(self.onselect_fasta_file)
        grid.addWidget(fasta_files_listbox, 3, 1)

        grid.addWidget(QLabel("Genes:"), 3, 2)
        genes_listbox = QListWidget()
        genes_listbox.setObjectName("genes_listbox")
        grid.addWidget(genes_listbox, 3, 3)


        rates_groupbox = QGroupBox()
        rates_groupbox_grid = QGridLayout()
        rates_groupbox.setLayout(rates_groupbox_grid)
        grid.addWidget(rates_groupbox, 3, 4)
        rates_groupbox_grid.addWidget(QLabel("Initiation rate: "), 0, 0)
        init_rate_spinbox = QDoubleSpinBox()
        init_rate_spinbox.setObjectName("init_rate_spinbox")
        init_rate_spinbox.setRange(0.01, 100.00)
        init_rate_spinbox.setSingleStep(0.01)
        rates_groupbox_grid.addWidget(init_rate_spinbox, 0, 1)

        rates_groupbox_grid.addWidget(QLabel("Termination rate: "), 1, 0)
        term_rate_spinbox = QDoubleSpinBox()
        term_rate_spinbox.setObjectName("term_rate_spinbox")
        term_rate_spinbox.setRange(0.01, 100.00)
        term_rate_spinbox.setSingleStep(0.01)
        rates_groupbox_grid.addWidget(term_rate_spinbox, 1, 1)

        rates_groupbox_grid.addWidget(QLabel("Transcript copy number: "), 2, 0)
        gene_copy_number_spinbox = QSpinBox()
        genes_listbox.setSelectionMode(2)
        gene_copy_number_spinbox.setObjectName("gene_copy_number_spinbox")
        gene_copy_number_spinbox.setRange(1, 1000)
        gene_copy_number_spinbox.setSingleStep(1)
        rates_groupbox_grid.addWidget(gene_copy_number_spinbox, 2, 1)

        add_gene_button = QPushButton("Add gene")
        add_gene_button.clicked.connect(self.add_simulation_entry)
        grid.addWidget(add_gene_button, 4, 1, 1, 4)

        grid.addWidget(QLabel("Added simulations: "), 5, 0)
        added_simulations_listbox = QTableWidget()
        added_simulations_listbox.setObjectName("added_simulations_listbox")
        added_simulations_listbox.setColumnCount(5)
        added_simulations_listbox.setShowGrid(False)
        added_simulations_listbox.setHorizontalHeaderLabels(['Fasta file', 'Gene',
                                                             'Initiation\nRate',
                                                             'Termination\nRate',
                                                             'Transcript\nCopy Number'])
        added_simulations_listbox.resizeColumnsToContents()
        added_simulations_listbox.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        added_simulations_listbox.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        grid.addWidget(added_simulations_listbox, 6, 0, 6, 6)

        remove_entry_button = QPushButton("Remove entry")
        remove_entry_button.clicked.connect(self.remove_simulation_entry)
        grid.addWidget(remove_entry_button, 14, 5)

        termination_condition_groupbox = QGroupBox()
        termination_condition_groupbox.setTitle("Stop condition")
        termination_condition_groupbox_grid = QGridLayout()
        termination_condition_groupbox.setLayout(termination_condition_groupbox_grid)
        grid.addWidget(termination_condition_groupbox, 15, 0, 2, 2)

        iteration_limit_radiobutton = QRadioButton("Iteration limit:")
        iteration_limit_radiobutton.setObjectName("iteration_limit_radiobutton")
        termination_condition_groupbox_grid.addWidget(iteration_limit_radiobutton, 0, 0)
        iteration_limit_spinbox = QSpinBox()
        iteration_limit_spinbox.setObjectName("iteration_limit_spinbox")
        iteration_limit_spinbox.setRange(0, int(10e10))
        iteration_limit_spinbox.valueChanged.connect(self.changed_iteration_limit)
        termination_condition_groupbox_grid.addWidget(iteration_limit_spinbox, 0, 1)

        time_limit_radiobutton = QRadioButton("Time limit:")
        time_limit_radiobutton.setObjectName("time_limit_radiobutton")
        termination_condition_groupbox_grid.addWidget(time_limit_radiobutton, 1, 0)
        time_limit_spinbox = QDoubleSpinBox()
        time_limit_spinbox.setObjectName("time_limit_spinbox")
        time_limit_spinbox.setRange(0, 10e10)
        time_limit_spinbox.valueChanged.connect(self.changed_time_limit_entry)
        termination_condition_groupbox_grid.addWidget(time_limit_spinbox, 1, 1)

        finished_ribosomes_limit_radiobutton = QRadioButton("Finished Ribosomes:")
        finished_ribosomes_limit_radiobutton.setObjectName("finished_ribosomes_limit_radiobutton")
        termination_condition_groupbox_grid.addWidget(finished_ribosomes_limit_radiobutton, 2, 0)
        finished_ribosomes_spinbox = QSpinBox()
        finished_ribosomes_spinbox.setObjectName("finished_ribosomes_spinbox")
        finished_ribosomes_spinbox.setRange(0, int(10e10))
        finished_ribosomes_spinbox.valueChanged.connect(self.changed_finished_ribosomes)
        termination_condition_groupbox_grid.addWidget(finished_ribosomes_spinbox, 2,1)

        history_groupbox = QGroupBox()
        history_groupbox.setFlat(True)
        history_groupbox_grid = QGridLayout()
        history_groupbox.setLayout(history_groupbox_grid)
        grid.addWidget(history_groupbox, 15, 3, 2, 1)

        history_groupbox_grid.addWidget(QLabel("Number of history entries: "), 0, 0)
        history_size_spinbox = QSpinBox()
        history_size_spinbox.setObjectName("history_size_spinbox")
        history_size_spinbox.setRange(1, int(10e15))
        history_size_spinbox.setValue(100000)
        history_groupbox_grid.addWidget(history_size_spinbox, 0, 1)

        generate_simulation_button = QPushButton("Generate Simulation file")
        generate_simulation_button.clicked.connect(self.generate_simulation_file)
        grid.addWidget(generate_simulation_button, 17, 1, 3, 4)
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
        initiation_rate_spinbox = self.window.findChild(QDoubleSpinBox, 'init_rate_spinbox')
        termination_rate_spinbox = self.window.findChild(QDoubleSpinBox, 'term_rate_spinbox')
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
        genes_list_indexes_to_be_updated = [i for i, x in enumerate(genes_list)
                                            if x in selected_genes]
        genes_list_to_be_added = [i for i, x in enumerate(selected_genes)
                                  if x not in genes_list]
        indexes_to_update = set(
            fasta_files_list_indexes).intersection(genes_list_indexes_to_be_updated)

        if len(indexes_to_update) > 0:
            # update gene copy number, initiation rate and termination rate
            for i in indexes_to_update:
                added_simulations_listbox.setItem(i,
                                                  2,
                                                  QTableWidgetItem(
                                                      str(initiation_rate_spinbox.value())))
                added_simulations_listbox.setItem(i,
                                                  3,
                                                  QTableWidgetItem(
                                                      str(termination_rate_spinbox.value())))
                added_simulations_listbox.setItem(i,
                                                  4,
                                                  QTableWidgetItem(
                                                      str(gene_copy_number_spinbox.value())))
        #process new entries
        for selected_gene_index in genes_list_to_be_added:
            # new entry
            indexes_to_update = added_simulations_listbox.rowCount()
            added_simulations_listbox.insertRow(indexes_to_update)
            added_simulations_listbox.setItem(indexes_to_update, 0,
                                              QTableWidgetItem(selected_fasta_file))
            added_simulations_listbox.setItem(indexes_to_update, 1,
                                              QTableWidgetItem(selected_genes[selected_gene_index]))
            added_simulations_listbox.setItem(indexes_to_update, 2,
                                              QTableWidgetItem(
                                                  str(initiation_rate_spinbox.value())))
            added_simulations_listbox.setItem(indexes_to_update, 3,
                                              QTableWidgetItem(
                                                  str(termination_rate_spinbox.value())))

            added_simulations_listbox.setItem(indexes_to_update, 4,
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

    def changed_iteration_limit(self):
        """
        select iteration limit as the user changed its value.
        """
        iteration_limit_radiobutton = self.window.findChild(QRadioButton, "iteration_limit_radiobutton")
        iteration_limit_radiobutton.setChecked(True)
    
    def changed_time_limit_entry(self):
        """
        select time limit as the user changed its value.
        """
        time_limit_radiobutton = self.window.findChild(QRadioButton, "time_limit_radiobutton")
        time_limit_radiobutton.setChecked(True)

    def changed_finished_ribosomes(self):
        """
        select finished ribosomes as the user changed its value.
        """
        finished_ribosomes_limit_radiobutton = self.window.findChild(QRadioButton, "finished_ribosomes_limit_radiobutton")
        finished_ribosomes_limit_radiobutton.setChecked(True)
    
    
    def generate_simulation_file(self):
        """
        This method should assemble the json and ask the user for a place to save the file.
        """
        return


class SimulationBuilder:
    """
    class to generate JSON config file
    """
    def __init__(self):
        """
        Create fields with default values.
        """
        self.concentration_file = ""
        self.pre_populate = True
        self.mRNA_entries = [] # will be a list of dictionaries
        self.terminated_ribosomes = 10000
        self.history_size = 10000
        self.results = {} # dictionary of dictionaries with lists

if __name__ == "__main__":
    Gui()
