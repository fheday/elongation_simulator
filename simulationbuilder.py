'''
Gui module to build and save a simulation set
'''

import  tkinter
from tkinter import filedialog
from Bio import SeqIO

class SimulationBuilder():
    """
    class to generate GUI.
    """
    def __init__(self):
        self.pre_populate = False
        self.concentrations_file_path = None
        self.genes_dict = {}

        self.window = tkinter.Tk()
        self.window.option_add("*Font", ("Times", 10))
        self.window.title("Simulation builder")

        self.concentrations_label_stringvar = tkinter.StringVar(self.window)
        self.concentrations_label_stringvar.set("concentration file: ")
        self.gene_copy_number_stringvar = tkinter.StringVar(self.window)
        self.gene_copy_number_stringvar.set("1")

        self.window.columnconfigure(1, weight=1)

        padx = 7
        pady = 5

        tkinter.Label(self.window, textvariable=self.concentrations_label_stringvar)\
            .grid(row=0, column=0, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        tkinter.Button(self.window, text="Open file", command=self.open_concentrations_file)\
            .grid(row=0, column=1, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)

        tkinter.Label(self.window, text="Initiation rate: ")\
            .grid(row=1, column=0, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        tkinter.Spinbox(self.window, from_=0.01, to=100.0, format="%.3f", increment=0.01, width=5)\
            .grid(row=1, column=1, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)

        tkinter.Label(self.window, text="Termination rate: ")\
            .grid(row=2, column=0, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        tkinter.Spinbox(self.window, from_=0.01, to=100.0, format="%.3f", increment=0.01, width=5)\
            .grid(row=2, column=1, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)

        tkinter.Label(self.window, text="Pre-populate: ")\
            .grid(row=3, column=0, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        tkinter.Checkbutton(self.window, variable=self.pre_populate)\
            .grid(row=3, column=1, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)

        tkinter.Label(self.window, text="Add FASTA file :")\
            .grid(row=4, column=0, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        tkinter.Button(self.window, text="Add file", command=self.open_fasta_file)\
            .grid(row=4, column=1, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)

        tkinter.Label(self.window, text="Loaded FASTA files:")\
            .grid(row=5, column=0, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        fasta_files_listbox = tkinter.Listbox(self.window, name="fasta_files_listbox")
        fasta_files_listbox.configure(exportselection=False)
        fasta_files_listbox.grid(row=5, column=1, padx=(padx),
                                 pady=pady, sticky=tkinter.E + tkinter.W)
        fasta_files_listbox.bind('<<ListboxSelect>>', self.onselect_fasta_file)
        tkinter.Label(self.window, text="Genes:").\
            grid(row=5, column=2, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        genes_listbox = tkinter.Listbox(self.window, name="genes_listbox")
        genes_listbox.configure(exportselection=False)
        genes_listbox.grid(row=5, column=3, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        tkinter.Label(self.window, text="Gene copy number: ")\
            .grid(row=5, column=4, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        tkinter.Spinbox(self.window, from_=1, to=1000, increment=1,
                        textvariable=self.gene_copy_number_stringvar, width=5)\
                    .grid(row=5, column=5, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)

        tkinter.Button(self.window, text="Add gene", command=self.add_simulation_entry)\
            .grid(row=6, column=1, columnspan=4, padx=(padx), pady=pady,
                  sticky=tkinter.E + tkinter.W)

        tkinter.Label(self.window, text="Added simulations: ")\
            .grid(row=7, column=0, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        tkinter.Listbox(self.window, name="added_simulations_listbox")\
            .grid(row=8, column=0, columnspan=6, padx=(padx), pady=pady,
                  sticky=tkinter.E + tkinter.W)

        tkinter.Button(self.window, text="Remove entry", command=self.remove_simulation_entry)\
            .grid(row=9, column=5, padx=(padx), pady=pady, sticky=tkinter.E + tkinter.W)
        self.show()

    def show(self):
        """
        bring GUI window to front and show it.
        """
        self.window.lift()
        self.window.mainloop()

    def open_concentrations_file(self):
        """
        open dialog to select concentration file.
        """
        concentrations_file_chooser = filedialog.askopenfilename(title="Select file",
                                                                 filetypes=(\
                                                                     ("csv files", "*.csv"),))
        if concentrations_file_chooser is None:
            return
        self.concentrations_file_path = concentrations_file_chooser
        self.concentrations_label_stringvar.set(
            "concentration file: " + str(concentrations_file_chooser))

    def open_fasta_file(self):
        """
        Open dialog to select FASTA file with genes to simulate.
        Also populate gene box.
        """
        fasta_file_chooser = filedialog.\
            askopenfilename(title="Select file", filetypes=(("txt files", "*.txt"),
                                                            ("FASTA file", "*.fasta")))
        if fasta_file_chooser is None or fasta_file_chooser in self.genes_dict.keys():
            return
        self.genes_dict[fasta_file_chooser] = \
            [record.id for record in SeqIO.parse(fasta_file_chooser, "fasta")]
        #update fasta files listbox
        fasta_files_listbox = self.window.children['fasta_files_listbox']
        if fasta_files_listbox.size() > 0:
            fasta_files_listbox.delete(0, tkinter.END)
        fasta_files_listbox.insert(tkinter.END, *self.genes_dict.keys())
        fasta_files_listbox.activate(tkinter.END)
        #update genes listbox
        genes_listbox = self.window.children["genes_listbox"]
        if genes_listbox.size() > 0:
            genes_listbox.delete(0, tkinter.END)
        genes_listbox.insert(tkinter.END, *self.genes_dict[fasta_file_chooser])
        genes_listbox.activate(tkinter.END)

    def onselect_fasta_file(self, evt):
        """
        Update gene box with the genes of the selected fasta file.
        """
        event_widget = evt.widget
        if len(event_widget.curselection()) == 0:
            return
        index = int(event_widget.curselection()[0])
        value = event_widget.get(index)
        #update genes listbox
        genes_listbox = self.window.children["genes_listbox"]
        if genes_listbox.size() > 0:
            genes_listbox.delete(0, tkinter.END)
        genes_listbox.insert(tkinter.END, *self.genes_dict[value])
        genes_listbox.activate(tkinter.END)

    def add_simulation_entry(self):
        """
        Use the selected fasta file, gene and gene copy number to add a new simulation entry.
        If fasta file AND gene is already added, update that entry with new copy number.
        """
        # we need: a gene and copy number.
        genes_listbox = self.window.children["genes_listbox"]
        fasta_files_listbox = self.window.children["fasta_files_listbox"]
        added_simulations_listbox = self.window.children["added_simulations_listbox"]
        if len(genes_listbox.curselection()) != 1 or len(fasta_files_listbox.curselection()) != 1:
            return
        selected_fasta_file = fasta_files_listbox.get(fasta_files_listbox.curselection()[0])
        selected_gene = genes_listbox.get(genes_listbox.curselection()[0])
        list_of_already_selected_entries = added_simulations_listbox.get(0, tkinter.END)
        update_entry = False
        index = 0
        for i, entry in enumerate(list_of_already_selected_entries):
            if selected_fasta_file + "    " + selected_gene in entry:
                update_entry = True
                index = i
                break
        if update_entry:
            #delete old entry
            added_simulations_listbox.delete(index)
            # add updated entry in same place
            added_simulations_listbox.insert(
                index, selected_fasta_file + "    " +\
                    selected_gene + "    " + self.gene_copy_number_stringvar.get())
        else:
            # new entry
            added_simulations_listbox.insert(tkinter.END,
                                             selected_fasta_file + "    " +\
                                             selected_gene + "    " +\
                                             self.gene_copy_number_stringvar.get())
        return

    def remove_simulation_entry(self):
        """
        Remove selected simulation entry.
        """
        added_simulations_listbox = self.window.children["added_simulations_listbox"]
        if len(added_simulations_listbox.curselection()) != 1:
            return
        added_simulations_listbox.delete(added_simulations_listbox.curselection()[0])
        return

if __name__ == "__main__":
    SimulationBuilder()
