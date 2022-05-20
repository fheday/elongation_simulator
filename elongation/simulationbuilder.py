'''
Gui module to build and save a simulation set
'''

import os
from PyQt5.QtWidgets import QApplication, QWidget, QGridLayout, QLabel, QPushButton,\
    QFileDialog, QSpinBox, QDoubleSpinBox, QCheckBox, QListWidget, QTableWidget,\
    QHeaderView, QTableWidgetItem, QGroupBox, QRadioButton
from PyQt5 import QtCore
from Bio import SeqIO
import pandas as pd
import json

class Gui():
    """
    class to generate GUI.
    """
    def __init__(self):
        self.genes_dict = {}
        self.stop_condition = ()

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
        add_fasta_file_button = QPushButton("Add FASTA file")
        add_fasta_file_button.setToolTip(\
            "Reads FASTA file and make its genes available for simulation")
        add_fasta_file_button.clicked.connect(self.open_fasta_file)
        grid.addWidget(add_fasta_file_button, 2, 1)

        load_batch_file_button = QPushButton("Load batch CSV file")
        load_batch_file_button.setToolTip("load a CSV file with a list of Genes, copy number, initiation and termination rates.")
        load_batch_file_button.clicked.connect(self.load_batch_file)
        grid.addWidget(load_batch_file_button, 2, 3)

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
        gene_copy_number_spinbox = QDoubleSpinBox()
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

        steady_state_checkbox = QCheckBox("After steady state: ")
        steady_state_checkbox.setObjectName("steady_state_checkbox")
        termination_condition_groupbox_grid.addWidget(steady_state_checkbox, 0, 0)
        steady_state_checkbox.clicked.connect(self.changed_steady_state_checkbox)
        # termination_condition_groupbox_grid.addWidget(steady_state_checkbox, 0, 1)

        iteration_limit_radiobutton.setObjectName("iteration_limit_radiobutton")
        termination_condition_groupbox_grid.addWidget(iteration_limit_radiobutton, 1, 0)
        iteration_limit_spinbox = QSpinBox()
        iteration_limit_spinbox.setObjectName("iteration_limit_spinbox")
        iteration_limit_spinbox.setRange(0, 2147483647)
        iteration_limit_spinbox.valueChanged.connect(self.changed_iteration_limit)
        termination_condition_groupbox_grid.addWidget(iteration_limit_spinbox, 1, 1)

        time_limit_radiobutton = QRadioButton("Time limit:")
        time_limit_radiobutton.setObjectName("time_limit_radiobutton")
        time_limit_radiobutton.setChecked(True) # have a checked option for starters.
        termination_condition_groupbox_grid.addWidget(time_limit_radiobutton, 2, 0)
        time_limit_spinbox = QDoubleSpinBox()
        time_limit_spinbox.setObjectName("time_limit_spinbox")
        time_limit_spinbox.setRange(0, 10e10)
        time_limit_spinbox.setValue(60) # set a default value
        self.stop_condition = ("time", float(time_limit_spinbox.value()))
        time_limit_spinbox.valueChanged.connect(self.changed_time_limit_entry)
        termination_condition_groupbox_grid.addWidget(time_limit_spinbox, 2, 1)

        finished_ribosomes_limit_radiobutton = QRadioButton("Terminated Ribosomes:")
        finished_ribosomes_limit_radiobutton.setObjectName("finished_ribosomes_limit_radiobutton")
        termination_condition_groupbox_grid.addWidget(finished_ribosomes_limit_radiobutton, 3, 0)
        finished_ribosomes_spinbox = QSpinBox()
        finished_ribosomes_spinbox.setObjectName("finished_ribosomes_spinbox")
        finished_ribosomes_spinbox.setRange(0, 2147483647)
        finished_ribosomes_spinbox.valueChanged.connect(self.changed_finished_ribosomes)
        termination_condition_groupbox_grid.addWidget(finished_ribosomes_spinbox, 3,1)

        history_groupbox = QGroupBox()
        history_groupbox.setFlat(True)
        history_groupbox_grid = QGridLayout()
        history_groupbox.setLayout(history_groupbox_grid)
        grid.addWidget(history_groupbox, 15, 3, 2, 1)

        history_groupbox_grid.addWidget(QLabel("Number of history entries: "), 0, 0)
        history_size_spinbox = QSpinBox()
        history_size_spinbox.setObjectName("history_size_spinbox")
        history_size_spinbox.setRange(1, 2147483647)
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

    def changed_steady_state_checkbox(self, state):
        steady_state_checkbox = self.window.findChild(QCheckBox, 'steady_state_checkbox')
        iteration_limit_radiobutton = self.window.findChild(QRadioButton, "iteration_limit_radiobutton")
        iteration_limit_spinbox = self.window.findChild(QSpinBox, "iteration_limit_spinbox")
        time_limit_radiobutton = self.window.findChild(QRadioButton, "time_limit_radiobutton")
        if iteration_limit_radiobutton.isChecked():
            time_limit_radiobutton.setChecked(True)
        iteration_limit_radiobutton.setChecked(False)
        iteration_limit_radiobutton.setEnabled(not steady_state_checkbox.checkState())
        iteration_limit_spinbox.setEnabled(not steady_state_checkbox.checkState())
        # update termination condition
        finished_ribosomes_limit_radiobutton = self.window.findChild(QRadioButton, "finished_ribosomes_limit_radiobutton")
        if iteration_limit_radiobutton.isChecked():
            self.changed_iteration_limit()
        elif time_limit_radiobutton.isChecked():
            self.changed_time_limit_entry()
        else:
            self.changed_finished_ribosomes()
            
        return

    def open_fasta_file(self):
        """
        Open dialog to select FASTA file with genes to simulate.
        Also populate gene box.
        """
        fasta_file_chooser, _ = QFileDialog\
            .getOpenFileName(self.window, "Select file", os.getcwd(),
                             "All compatible files (*.txt *.fasta *.fna *.fa);;txt File (*.txt );;FASTA file (*.fasta);;FNA file (*.fna *.fa)")
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
    
    def load_batch_file(self):
        """
        Open dialog to select CSV file with ORF, transcript copy number, initiation and termination rates.
        """
        batch_file_path, _ = QFileDialog.getOpenFileName(self.window, "Select file",
                                                         os.getcwd(), "CSV File (*.csv )")
        added_simulations_listbox = self.window.findChild(QTableWidget, 'added_simulations_listbox')
        if batch_file_path == '':
            return
        df = pd.read_csv(batch_file_path)
        if not set(['ORF', 'RNA_copies', 'Ini_rate', 'Term_rate']).issubset(df.columns):
            return
        for index, row in df.iterrows():
            for key in self.genes_dict.keys():
                index_to_update = added_simulations_listbox.rowCount()
                if row.loc["ORF"] in self.genes_dict[key]:
                    items = added_simulations_listbox.findItems(row.loc["ORF"], QtCore.Qt.MatchExactly)
                    if items:
                        index_to_update = items[0].row() # go to the row to edit entry.
                    else:
                        added_simulations_listbox.insertRow(index_to_update) # create new entry
                added_simulations_listbox.setItem(index_to_update, 0,
                                                QTableWidgetItem(key))
                added_simulations_listbox.setItem(index_to_update, 1,
                                                QTableWidgetItem(row.loc["ORF"]))
                added_simulations_listbox.setItem(index_to_update, 2,
                                                QTableWidgetItem(
                                                    str(row.loc["Ini_rate"])))
                added_simulations_listbox.setItem(index_to_update, 3,
                                                QTableWidgetItem(
                                                    str(row.loc["Term_rate"])))
                added_simulations_listbox.setItem(index_to_update, 4,
                                                QTableWidgetItem(
                                                    str(row.loc["RNA_copies"])))



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

    def get_entries_to_add_and_update(self, selected_fasta_file, selected_genes):
        added_simulations_listbox = self.window.findChild(QTableWidget, 'added_simulations_listbox')
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
        return genes_list_to_be_added, indexes_to_update

    def add_simulation_entry(self):
        """
        Use the selected fasta file, gene and gene copy number to add a new simulation entry.
        If fasta file AND gene is already added, update that entry with new copy number.
        """
        # we need: a gene and copy number.
        genes_listbox = self.window.findChild(QListWidget, 'genes_listbox')
        fasta_files_listbox = self.window.findChild(QListWidget, 'fasta_files_listbox')
        added_simulations_listbox = self.window.findChild(QTableWidget, 'added_simulations_listbox')
        gene_copy_number_spinbox = self.window.findChild(QDoubleSpinBox, 'gene_copy_number_spinbox')
        initiation_rate_spinbox = self.window.findChild(QDoubleSpinBox, 'init_rate_spinbox')
        termination_rate_spinbox = self.window.findChild(QDoubleSpinBox, 'term_rate_spinbox')
        if len(genes_listbox.selectedItems()) == 0 or len(fasta_files_listbox.selectedItems()) != 1:
            return
        selected_fasta_file = fasta_files_listbox.currentItem().text()
        selected_genes = [item.text() for item in genes_listbox.selectedItems()]
        genes_list_to_be_added, indexes_to_update = self.get_entries_to_add_and_update(selected_fasta_file, selected_genes)

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
        iteration_limit_spinbox = self.window.findChild(QSpinBox, "iteration_limit_spinbox")
        iteration_limit_radiobutton.setChecked(True)
        self.stop_condition = ("iteration", int(iteration_limit_spinbox.value()))
    
    def changed_time_limit_entry(self):
        """
        select time limit as the user changed its value.
        """
        time_limit_radiobutton = self.window.findChild(QRadioButton, "time_limit_radiobutton")
        time_limit_spinbox = self.window.findChild(QDoubleSpinBox, "time_limit_spinbox")
        steady_state_checkbox = self.window.findChild(QCheckBox, 'steady_state_checkbox')
        time_limit_radiobutton.setChecked(True)
        if steady_state_checkbox.checkState():
            self.stop_condition = ("steady_state_time", float(time_limit_spinbox.value()))
        else:
            self.stop_condition = ("time", float(time_limit_spinbox.value()))

    def changed_finished_ribosomes(self):
        """
        select finished ribosomes as the user changed its value.
        """
        finished_ribosomes_limit_radiobutton = self.window.findChild(QRadioButton, "finished_ribosomes_limit_radiobutton")
        finished_ribosomes_spinbox = self.window.findChild(QSpinBox, "finished_ribosomes_spinbox")
        steady_state_checkbox = self.window.findChild(QCheckBox, 'steady_state_checkbox')
        finished_ribosomes_limit_radiobutton.setChecked(True)
        if steady_state_checkbox.checkState():
            self.stop_condition = ("steady_state_ribosomes", int(finished_ribosomes_spinbox.value()))
        else:
            self.stop_condition = ("ribosomes", int(finished_ribosomes_spinbox.value()))
    
    
    def generate_simulation_file(self):
        """
        This method should assemble the json and ask the user for a place to save the file.
        """
        file_name = QFileDialog.getSaveFileName(self.window, 'Save simulation configuration', os.getcwd(), 'Simulation Files (*.json) ;;All files(*)', options=QFileDialog.DontUseNativeDialog)
        if file_name == "":
            return
        if file_name[0] and file_name[1].split('.')[1] in file_name[0].split(".")[-1]:
            file_name = file_name[0] # file name already has extension.
        else:
            file_name = file_name[0] + '.' + file_name[1].split('.')[1][:-1]
        selected_concentration_label = self.window.findChild(QLabel, "selected_concentration_label")
        pre_populate_check_box = self.window.findChild(QCheckBox, "pre_populate_check_box")
        history_size_spinbox = self.window.findChild(QSpinBox, "history_size_spinbox")
        added_simulations_table = self.window.findChild(QTableWidget, 'added_simulations_listbox')
        
        sb = SimulationBuilder()
        sb.set_concentration_file(selected_concentration_label.text())
        sb.set_pre_populate(pre_populate_check_box.checkState())
        sb.set_history_size(int(history_size_spinbox.value()))
        sb.set_halt_condition(*self.stop_condition)
        for index in range(added_simulations_table.rowCount()):
            fasta_file = added_simulations_table.item(index, 0).text()
            gene = added_simulations_table.item(index, 1).text()
            initiation_rate = float(added_simulations_table.item(index, 2).text())
            termination_rate = float(added_simulations_table.item(index, 3).text())
            transcript_copy_number = float(added_simulations_table.item(index, 4).text())
            sb.add_mRNA_entry(fasta_file, gene, initiation_rate, termination_rate, transcript_copy_number)
        sb.save_simulation(file_name)

        return


class SimulationBuilder:
    """
    class to generate JSON config file
    """
    def __init__(self):
        """
        Create fields with default values.
        """
        self.data = {
            "concentration_file":"",
            "pre_populate":True,
            "mRNA_entries":[],
            "history_size":100000,
        }

    def set_concentration_file(self, file_path):
        self.data["concentration_file"] = file_path
    
    def set_pre_populate(self, pre_populate):
        self.data["pre_populate"] = pre_populate

    def set_history_size(self, history_size):
        if history_size > 0:
            self.data["history_size"] = history_size
    
    def add_mRNA_entry(self, fasta_file, gene, initiation_rate, termination_rate, transcript_copy_number):
        if initiation_rate < 0 or termination_rate < 0 or transcript_copy_number <= 0:
            return
        entry = {
            "fasta_file":fasta_file,
            "gene":gene,
            "initiation_rate":initiation_rate,
            "termination_rate":termination_rate,
            "transcript_copy_number":transcript_copy_number
        }
        if entry not in self.data["mRNA_entries"]:
            self.data["mRNA_entries"].append(entry)

    def set_halt_condition(self, condition, value):
        if value < 0:
            return
        if condition.lower() == "iteration":
            self.data["iteration_limit"] = int(value)
        elif condition.lower() == "time":
            self.data["time_limit"] = value
        elif condition.lower() == "ribosomes":
            self.data["finished_ribosomes"] = int(value)
        elif condition.lower() == "steady_state_time":
            self.data["steady_state_time"] = value
        elif condition.lower() == "steady_state_ribosomes":
            self.data["steady_state_ribosomes"] = value
    
    def save_simulation(self, file_path):
        # validation
        if len(set(["iteration_limit", "time_limit", "finished_ribosomes", "steady_state_time", "steady_state_ribosomes"]).intersection(self.data.keys())) == 0:
            return
        if len(self.data["mRNA_entries"]) == 0:
            return
        # convert absolute paths to relative path
        base_path = os.path.split(file_path)[0] + os.sep
        self.data["concentration_file"] = os.path.relpath(self.data["concentration_file"], start=base_path)
        for mRNA_entry in self.data["mRNA_entries"]:
            mRNA_entry["fasta_file"] = os.path.relpath(mRNA_entry["fasta_file"], start=base_path)
        # valid simulation. save.
        with open(file_path, 'w') as outfile:
            json.dump(self.data, outfile, indent=4)


if __name__ == "__main__":
    Gui()
