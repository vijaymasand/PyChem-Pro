# PyChem-Pro — pdbx/writer/PdbxWriter.py
#
# PDBx/mmCIF serialiser — writes DataContainer objects to CIF format.
# Migrated from Python 2 to Python 3 by the PyChem-Pro project (2026-06-21).
#   Changes: fixed integer division (rowCount / rowPartition → //),
#            added UTF-8 encoding guidance, removed Python 2 print statements.
#
"""
PDBx/mmCIF writer — serialises DataContainer objects to mmCIF text.

Usage::

    from pdbx.reader.PdbxContainers import DataContainer, DataCategory
    from pdbx.writer.PdbxWriter import PdbxWriter

    with open("out.cif", "w", encoding="utf-8") as fh:
        writer = PdbxWriter(fh)
        writer.write(container_list)
"""

__version__ = "3.0.0"

import sys
from pdbx.reader.PdbxContainers import DataContainer, DataCategory, DefinitionContainer



class PdbxError(Exception):
    """Class for catching general PDBx errors."""
    pass


class PdbxWriter(object):
    """Write PDBx/mmCIF data files or dictionaries.

    Accepts a list of DataContainer / DefinitionContainer objects and
    serialises them to the output file handle in compliant mmCIF syntax.

    Args:
        ofh: Output file handle (default: sys.stdout).
              Open with ``encoding='utf-8'`` for non-ASCII data.
    """

    def __init__(self, ofh=sys.stdout):
        self.__ofh = ofh
        self.__containerList = []
        self.__MAXIMUM_LINE_LENGTH = 2048
        self.__SPACING = 2
        self.__INDENT_DEFINITION = 3
        self.__indentSpace = " " * self.__INDENT_DEFINITION
        self.__doDefinitionIndent = False
        # Maximum number of rows sampled for value-length/format scanning.
        self.__rowPartition = None

    def setRowPartition(self, numRows):
        """Set the maximum number of rows checked for value length and format.

        For very large loop_ tables this provides a significant speed-up.
        """
        self.__rowPartition = numRows

    def write(self, containerList):
        """Write all containers to the output file handle."""
        self.__containerList = containerList
        for container in self.__containerList:
            self.writeContainer(container)

    def writeContainer(self, container):
        """Write a single data or definition container."""
        indS = " " * self.__INDENT_DEFINITION

        if isinstance(container, DefinitionContainer):
            self.__write("save_%s\n" % container.getName())
            self.__doDefinitionIndent = True
            self.__write(indS + "#\n")
        elif isinstance(container, DataContainer):
            if container.getGlobal():
                self.__write("global_\n")
                self.__doDefinitionIndent = False
                self.__write("\n")
            else:
                self.__write("data_%s\n" % container.getName())
                self.__doDefinitionIndent = False
                self.__write("#\n")

        for nm in container.getObjNameList():
            obj = container.getObj(nm)
            objL = obj.getRowList()

            if len(objL) == 0:
                continue
            elif len(objL) == 1:
                self.__writeItemValueFormat(obj)
            elif len(objL) > 1 and len(obj.getAttributeList()) > 0:
                self.__writeTableFormat(obj)
            else:
                raise PdbxError(
                    "Unexpected state for category %s: %d rows, %d attributes" %
                    (nm, len(objL), len(obj.getAttributeList())))

            if self.__doDefinitionIndent:
                self.__write(indS + "#")
            else:
                self.__write("#")

        if isinstance(container, DefinitionContainer):
            self.__write("\nsave_\n")
        self.__write("#\n")

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def __write(self, st):
        self.__ofh.write(st)

    def __writeItemValueFormat(self, myCategory):
        """Write key-value pairs (single-row category)."""
        attributeNameLengthMax = 0
        for attributeName in myCategory.getAttributeList():
            attributeNameLengthMax = max(attributeNameLengthMax, len(attributeName))
        itemNameLengthMax = (self.__SPACING + len(myCategory.getName())
                             + attributeNameLengthMax + 2)

        lineList = ["#\n"]
        for attributeName, _iPos in myCategory.getAttributeListWithOrder():
            if self.__doDefinitionIndent:
                lineList.append(self.__indentSpace)
            itemName = "_%s.%s" % (myCategory.getName(), attributeName)
            lineList.append(itemName.ljust(itemNameLengthMax))
            lineList.append(myCategory.getValueFormatted(attributeName, 0))
            lineList.append("\n")

        self.__write("".join(lineList))

    def __writeTableFormat(self, myCategory):
        """Write a loop_ block (multi-row category)."""
        # ── loop_ header ──────────────────────────────────────────────
        lineList = ["#\n"]
        if self.__doDefinitionIndent:
            lineList.append(self.__indentSpace)
        lineList.append("loop_")
        for attributeName in myCategory.getAttributeList():
            lineList.append("\n")
            if self.__doDefinitionIndent:
                lineList.append(self.__indentSpace)
            lineList.append("_%s.%s" % (myCategory.getName(), attributeName))
        self.__write("".join(lineList))

        # ── Determine format types (optionally on a sample) ───────────
        if self.__rowPartition is not None:
            numSteps = max(1, myCategory.getRowCount() // self.__rowPartition)
        else:
            numSteps = 1

        formatTypeList, _dataTypeList = myCategory.getFormatTypeList(steps=numSteps)
        maxLengthList = myCategory.getAttributeValueMaxLengthList(steps=numSteps)
        spacing = " " * self.__SPACING

        # ── Loop data rows ────────────────────────────────────────────
        for iRow in range(myCategory.getRowCount()):
            lineList = ["\n"]
            if self.__doDefinitionIndent:
                lineList.append(self.__indentSpace + "  ")

            for iAt in range(myCategory.getAttributeCount()):
                formatType = formatTypeList[iAt]
                maxLength = maxLengthList[iAt]

                if formatType in ('FT_UNQUOTED_STRING', 'FT_NULL_VALUE'):
                    val = myCategory.getValueFormattedByIndex(iAt, iRow)
                    lineList.append(val.ljust(maxLength))

                elif formatType == 'FT_NUMBER':
                    val = myCategory.getValueFormattedByIndex(iAt, iRow)
                    lineList.append(val.rjust(maxLength))

                elif formatType == 'FT_QUOTED_STRING':
                    val = myCategory.getValueFormattedByIndex(iAt, iRow)
                    lineList.append(val.ljust(maxLength + 2))

                elif formatType == 'FT_MULTI_LINE_STRING':
                    val = myCategory.getValueFormattedByIndex(iAt, iRow)
                    lineList.append(val)

                lineList.append(spacing)

            self.__write("".join(lineList))

        self.__write("\n")
