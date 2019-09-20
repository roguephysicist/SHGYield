# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'InputParameter.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_bandgap(object):
    def setupUi(self, bandgap):
        bandgap.setObjectName("bandgap")
        bandgap.resize(515, 89)
        bandgap.setMinimumSize(QtCore.QSize(515, 74))
        bandgap.setMaximumSize(QtCore.QSize(515, 16777215))
        self.gridLayout = QtWidgets.QGridLayout(bandgap)
        self.gridLayout.setObjectName("gridLayout")
        self.doubleSpinBox = QtWidgets.QDoubleSpinBox(bandgap)
        self.doubleSpinBox.setObjectName("doubleSpinBox")
        self.gridLayout.addWidget(self.doubleSpinBox, 2, 4, 1, 1)
        self.label_5 = QtWidgets.QLabel(bandgap)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 0, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(bandgap)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 1, 1, 1)
        self.line = QtWidgets.QFrame(bandgap)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayout.addWidget(self.line, 1, 0, 1, 5)
        self.label_3 = QtWidgets.QLabel(bandgap)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 2, 2, 1, 1)
        self.label = QtWidgets.QLabel(bandgap)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 2, 0, 1, 1)
        self.label_7 = QtWidgets.QLabel(bandgap)
        self.label_7.setObjectName("label_7")
        self.gridLayout.addWidget(self.label_7, 0, 2, 1, 1)
        self.label_6 = QtWidgets.QLabel(bandgap)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 0, 1, 1, 1)
        self.label_4 = QtWidgets.QLabel(bandgap)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 0, 4, 1, 1)

        self.retranslateUi(bandgap)
        QtCore.QMetaObject.connectSlotsByName(bandgap)

    def retranslateUi(self, bandgap):
        _translate = QtCore.QCoreApplication.translate
        bandgap.setWindowTitle(_translate("bandgap", "Band Gap (Automated Scissors)"))
        self.doubleSpinBox.setSuffix(_translate("bandgap", " eV"))
        self.label_5.setText(_translate("bandgap", "Case"))
        self.label_2.setText(_translate("bandgap", "TextLabel"))
        self.label_3.setText(_translate("bandgap", "TextLabel"))
        self.label.setText(_translate("bandgap", "case"))
        self.label_7.setText(_translate("bandgap", "Scissors ℏΔ (eV)"))
        self.label_6.setText(_translate("bandgap", "Uncorrected (eV)"))
        self.label_4.setText(_translate("bandgap", "Target (eV)"))

