import sys, os, json
from dataclasses import dataclass
from typing import Iterable, Callable
from PyQt5.QtWidgets import (
    QApplication,
    QGridLayout,
    QTreeWidget,
    QTreeWidgetItem,
    QWidget,
    QMainWindow,
    QDialog,
    QLabel,
    QDialogButtonBox,
    QTextBrowser,
    QVBoxLayout,
    QLineEdit,
    QMessageBox,
    QPushButton
)
from PyQt5 import QtCore, QtGui

@dataclass
class Pairing_relationship:
    codon: str
    anticodons: list[str]

    @classmethod
    def load_dict_one_position(cls, json_file_name: str, pairing_type: str) -> list['Pairing_relationship']:
        result = []
        with open(json_file_name) as json_file:
            data = json.load(json_file)
            for codon in data[pairing_type]:
                result.append(cls(codon, data[pairing_type][codon]))
            return result
    
    @classmethod
    def convert_one_position_from_dict(cls, data: dict, pairing_type: str) -> list['Pairing_relationship']:
        result = []
        for codon in data[pairing_type]:
            result.append(cls(codon, data[pairing_type][codon]))
        return result


    @classmethod
    def load_all_dict(cls, json_file_name: str):
        result = {}
        if os.path.isfile(json_file_name) and os.path.exists(json_file_name):
            # this is an existing file. Edit:
            with open(json_file_name) as json_file:
                data = json.load(json_file)
                for pairing in data.keys():
                    if pairing == "Pairing Rules":
                        result[pairing] = data[pairing]
                    else:
                        result[pairing] = cls.load_dict_one_position(json_file_name, pairing)
        else:
            # need to create a new file.
            data = {"Watson-Crick": {"A": ["U", "&", "3", "1", "~", "N", "S", ")", "{", "V", "}", "P"],
                                     "C": ["G", "#", "W"], "G": ["C", "B"], "U": ["A"]},
                    "Wobble": {"A": ["A", "I", "M", "?"], "C": ["A", "U", "P", "I", "?", "Q"], 
                               "G": ["A", "U", "&", "3", "1", "~", "N", "S", ")", "{", "V", "P", "?", "M"],
                               "U": ["G", "#", "W", "U", "V","P", "I", "Q"]}}
            for pairing in data.keys():
                result[pairing] = cls.convert_one_position_from_dict(data, pairing)
        return result

@dataclass
class Classification_rules:
    classification_type: str
    rules: list[(str, str, str)]


class Window(QMainWindow):

    def __init__(self, pairings, display_file_name):
        super().__init__()
        self.setWindowTitle("Base pairing editor: " + display_file_name)
        widget = QWidget()
        layout = QGridLayout()

        tree = QTreeWidget()
        tree.setColumnCount(2)
        tree.setHeaderLabels(["Codon", "Anticodons"])
        tree.setColumnWidth(0, 130)

        root = QTreeWidgetItem(tree)
        root.setText(0, "Pairing type")

        self.watson_crick_pairing_widgetItem = QTreeWidgetItem(root)
        self.watson_crick_pairing_widgetItem.setText(0, "Watson-Crick")

        # fill Watson-Crick
        for pairing in pairings["Watson-Crick"]:
            node = QTreeWidgetItem((pairing.codon, ','.join(pairing.anticodons)))
            self.watson_crick_pairing_widgetItem.addChild(node)


        self.wobble_pairing_widgetItem = QTreeWidgetItem(root)
        self.wobble_pairing_widgetItem.setText(0, "Wobble")

        self.pairing_rules = QTreeWidgetItem(tree)
        self.pairing_rules.setText(0, "Pairing Rules")
        self.near_cognate = QTreeWidgetItem(self.pairing_rules)
        self.near_cognate.setText(0, "Near-Cognate")
        self.near_cognate_base_rules = QTreeWidgetItem(self.near_cognate)
        self.near_cognate_base_rules.setText(0, "base-level")
        self.near_cognate_base_rules.addChild(QTreeWidgetItem(("base-level", str(pairings["Pairing Rules"]["Near-Cognate"]["base-level"]))))

        
        # fill Wobble
        for pairing in pairings["Wobble"]:
            node = QTreeWidgetItem((pairing.codon, ','.join(pairing.anticodons)))
            self.wobble_pairing_widgetItem.addChild(node)

        tree.itemDoubleClicked.connect(self.edit_pairing)
        
        layout.addWidget(tree, 0, 0)

        self.save_button = QPushButton("Save")
        self.save_button.clicked.connect(self.save)
        layout.addWidget(self.save_button, 1,0)


        widget.setLayout(layout)
        self.setCentralWidget(widget)
        tree.expandAll()

        self.show()

    def identify_click(self, clicked_obj) -> (str, str):
        for child_number in range(self.watson_crick_pairing_widgetItem.childCount()):
            if self.watson_crick_pairing_widgetItem.child(child_number) == clicked_obj:
                return "Watson-Crick"
            if self.wobble_pairing_widgetItem.child(child_number) == clicked_obj:
                return "Wobble"
        return None
        
    def edit_pairing(self, clicked_obj):
        pairing_type = self.identify_click(clicked_obj)
        if pairing_type == None:
            return
        # open dialog for changing data.
        selected_pairing = Pairing_relationship(clicked_obj.data(0,0), clicked_obj.data(1,0))
        pairing_dialog = EditPairingDialog(self)
        pairing_dialog.set_pair(selected_pairing)
        
        pairing_dialog.exec_()
        if pairing_dialog.pairing_relationship.codon == selected_pairing.codon and\
           pairing_dialog.pairing_relationship.anticodons == selected_pairing.anticodons:
            # nothing changed.
            return
        else:
            # we need to update the tree.abs
            pairing_node = None
            if pairing_type == "Watson-Crick":
                pairing_node = self.watson_crick_pairing_widgetItem
            else:
                pairing_node = self.wobble_pairing_widgetItem
            for child_no in range(pairing_node.childCount()):
                item = pairing_node.child(child_no)
                if pairing_dialog.pairing_relationship.codon == item.data(0, 0):
                    # update
                    item.setData(0, 0, pairing_dialog.pairing_relationship.codon)
                    item.setData(1, 0, pairing_dialog.pairing_relationship.anticodons)

    def save(self):
        # recreate the dictionary
        data = {'Watson-Crick':{}, 'Wobble': {}}
        # Watson-Crick:
        for child_no in range(self.watson_crick_pairing_widgetItem.childCount()):
            item = self.watson_crick_pairing_widgetItem.child(child_no)
            data['Watson-Crick'][item.data(0, 0)] = [k for k in item.data(1,0).replace(' ', '').split(',')]
        
        # Wobble
        for child_no in range(self.wobble_pairing_widgetItem.childCount()):
            item = self.wobble_pairing_widgetItem.child(child_no)
            data['Wobble'][item.data(0, 0)] = [k for k in item.data(1,0).replace(' ', '').split(',')]
        # confirm save the data.
        qm = QMessageBox
        ret = qm.question(self,'', "Are you sure to save the file?", qm.Yes | qm.No)
        if ret == qm.Yes:
            # do save it
            with open(sys.argv[1], "w") as outfile:
                outfile.write(json.dumps(data, indent=4))

class EditPairingDialog(QDialog):
    def __init__(self, parent=None):
        super(EditPairingDialog, self).__init__(parent)
        self.setWindowTitle("Edit base pairing")
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.clicked_ok)
        self.buttonBox.rejected.connect(self.clicked_cancel)

        
        self.codon_textbox = QLineEdit(self)
        self.codon_textbox.setInputMask("A")

        self.anticodons_textbox = QLineEdit(self)
        self.anticodons_textbox.setFixedWidth(170)
        # self.validator = QtGui.QRegExpValidator(QtCore.QRegExp("+(\S)(,|$)"), self.anticodons_textbox)
        # self.anticodons_textbox.setValidator(self.validator)


        self.layout = QGridLayout(self)
        self.layout.addWidget(QLabel("Codon"), 0, 0)
        self.layout.addWidget(self.codon_textbox, 0, 1)
        self.layout.addWidget(QLabel("Anticodons"), 1, 0)
        self.layout.addWidget(self.anticodons_textbox, 1, 1)
        self.layout.addWidget(self.buttonBox, 2, 0)

    def clicked_ok(self):
        if all( len(k)==1 for k in self.anticodons_textbox.text().replace(' ', '').split(',')):
            # update the pairing
            self.pairing_relationship = Pairing_relationship(self.codon_textbox.text(), self.anticodons_textbox.text())
            self.close()
        else:
            QMessageBox.about(self, "Invalid input", "The input is invalid: Either there is a codon with a relationship or the list of anitcodons is not comma separated. Please review your input.")
        

    def clicked_cancel(self):
        self.close()

    def set_pair(self, pair):
        self.pairing_relationship = pair
        self.codon_textbox.setText(self.pairing_relationship.codon)
        self.anticodons_textbox.setText(self.pairing_relationship.anticodons)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    pairings = None
    file_name = sys.argv[1]
    if len(sys.argv) > 1:
        pairings = Pairing_relationship.load_all_dict(file_name)
    else:
        raise ValueError("Base pairing file not informed.")
    window = Window(pairings, os.path.basename(file_name))
    window.resize(int(window.width()*1.5), int(window.height()*2.5))
    sys.exit(app.exec_())