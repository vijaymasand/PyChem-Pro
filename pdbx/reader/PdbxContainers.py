# PyChem-Pro — pdbx/reader/PdbxContainers.py
#
# Container classes supporting the PDBx/mmCIF storage model.
# Migrated from Python 2 to Python 3 by the PyChem-Pro project (2026-06-21).
#   Changes: xrange → range, fixed setValue / invokeAttributeMethod /
#            getAttributeLengthMaximumList, removed Python 2 has_key() usage,
#            type-checking modernised (isinstance), silent bare-except cleaned.
#
"""
Container classes for the PDBx/mmCIF storage model.

Provides DataContainer, DataCategory, and supporting base classes that
represent data blocks and tabular categories in PDBx/mmCIF files.

A base container class supports common features of data and definition
containers.  PDBx data files are organised in sections called data blocks
which are mapped to DataContainer objects.  PDBx dictionaries contain
definition sections and data sections mapped to DefinitionContainer and
DataContainer objects respectively.

Data categories use '_categoryName.attributeName' item labels (analogous
to table.column in SQL, or class.attribute in OO design).

The DataCategory class provides the primary storage container for
instance data and definition metadata.
"""

__version__ = "3.0.0"

import re
import sys
import traceback



class CifName(object):
    """Class of utilities for CIF-style data names."""

    def __init__(self):
        pass

    @staticmethod
    def categoryPart(name):
        tname = ""
        if name.startswith("_"):
            tname = name[1:]
        else:
            tname = name
        i = tname.find(".")
        if i == -1:
            return tname
        else:
            return tname[:i]

    @staticmethod
    def attributePart(name):
        i = name.find(".")
        if i == -1:
            return None
        else:
            return name[i + 1:]


class ContainerBase(object):
    """Container base class for data and definition objects."""

    def __init__(self, name):
        # The enclosing scope of the data container (e.g. data_/save_)
        self.__name = name
        # List of category names within this container
        self.__objNameList = []
        # Dictionary of DataCategory objects keyed by category name
        self.__objCatalog = {}
        self.__type = None

    def getType(self):
        return self.__type

    def setType(self, type):
        self.__type = type

    def getName(self):
        return self.__name

    def setName(self, name):
        self.__name = name

    def exists(self, name):
        return name in self.__objCatalog

    def getObj(self, name):
        if name in self.__objCatalog:
            return self.__objCatalog[name]
        return None

    def getObjNameList(self):
        return self.__objNameList

    def append(self, obj):
        """Add the input object to the current object catalog.
        An existing object of the same name will be overwritten.
        """
        if obj.getName() is not None:
            if obj.getName() not in self.__objCatalog:
                self.__objNameList.append(obj.getName())
            self.__objCatalog[obj.getName()] = obj

    def replace(self, obj):
        """Replace an existing object with the input object."""
        if (obj.getName() is not None) and (obj.getName() in self.__objCatalog):
            self.__objCatalog[obj.getName()] = obj

    def printIt(self, fh=sys.stdout, type="brief"):
        fh.write("+ %s container: %30s contains %4d categories\n" %
                 (self.getType(), self.getName(), len(self.__objNameList)))
        for nm in self.__objNameList:
            fh.write("--------------------------------------------\n")
            fh.write("Data category: %s\n" % nm)
            if type == 'brief':
                self.__objCatalog[nm].printIt(fh)
            else:
                self.__objCatalog[nm].dumpIt(fh)

    def rename(self, curName, newName):
        """Change the name of an object in place."""
        try:
            i = self.__objNameList.index(curName)
            self.__objNameList[i] = newName
            self.__objCatalog[newName] = self.__objCatalog[curName]
            self.__objCatalog[newName].setName(newName)
            del self.__objCatalog[curName]
            return True
        except Exception:
            return False

    def remove(self, curName):
        """Remove object by name. Return True on success or False otherwise."""
        try:
            if curName in self.__objCatalog:
                del self.__objCatalog[curName]
                i = self.__objNameList.index(curName)
                del self.__objNameList[i]
                return True
        except Exception:
            pass
        return False


class DefinitionContainer(ContainerBase):
    def __init__(self, name):
        super(DefinitionContainer, self).__init__(name)
        self.setType('definition')

    def isCategory(self):
        return self.exists('category')

    def isAttribute(self):
        return self.exists('item')

    def printIt(self, fh=sys.stdout, type="brief"):
        fh.write("Definition container: %30s contains %4d categories\n" %
                 (self.getName(), len(self.getObjNameList())))
        if self.isCategory():
            fh.write("Definition type: category\n")
        elif self.isAttribute():
            fh.write("Definition type: item\n")
        else:
            fh.write("Definition type: undefined\n")
        for nm in self.getObjNameList():
            fh.write("--------------------------------------------\n")
            fh.write("Definition category: %s\n" % nm)
            if type == 'brief':
                self.getObj(nm).printIt(fh)
            else:
                self.getObj(nm).dumpIt(fh)


class DataContainer(ContainerBase):
    """Container class for DataCategory objects."""

    def __init__(self, name):
        super(DataContainer, self).__init__(name)
        self.setType('data')
        self.__globalFlag = False

    def invokeDataBlockMethod(self, type, method, db):
        self.__currentRow = 1
        exec(method.getInline())

    def setGlobal(self):
        self.__globalFlag = True

    def getGlobal(self):
        return self.__globalFlag


class DataCategoryBase(object):
    """Base object definition for a data category."""

    def __init__(self, name, attributeNameList=None, rowList=None):
        self._name = name
        self._rowList = rowList if rowList is not None else []
        self._attributeNameList = attributeNameList if attributeNameList is not None else []
        self._catalog = {}
        self._numAttributes = 0
        self.__setup()

    def __setup(self):
        self._numAttributes = len(self._attributeNameList)
        self._catalog = {}
        for attributeName in self._attributeNameList:
            self._catalog[attributeName.lower()] = attributeName

    def setRowList(self, rowList):
        self._rowList = rowList

    def setAttributeNameList(self, attributeNameList):
        self._attributeNameList = attributeNameList
        self.__setup()

    def setName(self, name):
        self._name = name

    def get(self):
        return (self._name, self._attributeNameList, self._rowList)


class DataCategory(DataCategoryBase):
    """Methods for creating, accessing, and formatting PDBx cif data categories."""

    def __init__(self, name, attributeNameList=None, rowList=None):
        super(DataCategory, self).__init__(name, attributeNameList, rowList)
        self.__lfh = sys.stdout
        self.__currentRowIndex = 0
        self.__currentAttribute = None
        self.__avoidEmbeddedQuoting = False

        # Regular expressions for quoting decisions
        self.__wsRe = re.compile(r"\s")
        self.__wsAndQuotesRe = re.compile(r"[\s'\"]")
        self.__nlRe = re.compile(r"[\n\r]")
        self.__sqRe = re.compile(r"[']")
        self.__sqWsRe = re.compile(r"('\s)|(\s')")
        self.__dqRe = re.compile(r'["]')
        self.__dqWsRe = re.compile(r'("\s)|(\s")')
        self.__intRe = re.compile(r'^[0-9]+$')
        self.__floatRe = re.compile(
            r'^-?(([0-9]+)[.]?|([0-9]*[.][0-9]+))([(][0-9]+[)])?([eE][+-]?[0-9]+)?$')

        self.__dataTypeList = [
            'DT_NULL_VALUE', 'DT_INTEGER', 'DT_FLOAT', 'DT_UNQUOTED_STRING',
            'DT_ITEM_NAME', 'DT_DOUBLE_QUOTED_STRING', 'DT_SINGLE_QUOTED_STRING',
            'DT_MULTI_LINE_STRING'
        ]
        self.__formatTypeList = [
            'FT_NULL_VALUE', 'FT_NUMBER', 'FT_NUMBER', 'FT_UNQUOTED_STRING',
            'FT_QUOTED_STRING', 'FT_QUOTED_STRING', 'FT_QUOTED_STRING',
            'FT_MULTI_LINE_STRING'
        ]

    def __getitem__(self, x):
        """Implements list-type functionality.

        x=integer  -> returns the row in category (normal list behavior)
        x=string   -> returns the value of attribute 'x' in first row.
        """
        if isinstance(x, int):
            return self._rowList[x]
        elif isinstance(x, str):
            try:
                ii = self.getAttributeIndex(x)
                return self._rowList[0][ii]
            except (IndexError, KeyError):
                raise KeyError
        raise TypeError(x)

    def getCurrentAttribute(self):
        return self.__currentAttribute

    def getRowIndex(self):
        return self.__currentRowIndex

    def getRowList(self):
        return self._rowList

    def getRowCount(self):
        return len(self._rowList)

    def getRow(self, index):
        try:
            return self._rowList[index]
        except Exception:
            return []

    def removeRow(self, index):
        try:
            if 0 <= index < len(self._rowList):
                del self._rowList[index]
                if self.__currentRowIndex >= len(self._rowList):
                    self.__currentRowIndex = max(0, len(self._rowList) - 1)
                return True
        except Exception:
            pass
        return False

    def getFullRow(self, index):
        """Return a full row based on the length of the attribute list."""
        try:
            if len(self._rowList[index]) < self._numAttributes:
                for ii in range(self._numAttributes - len(self._rowList[index])):
                    self._rowList[index].append('?')
            return self._rowList[index]
        except Exception:
            return ['?' for _ in range(self._numAttributes)]

    def getName(self):
        return self._name

    def getAttributeList(self):
        return self._attributeNameList

    def getAttributeCount(self):
        return len(self._attributeNameList)

    def getAttributeListWithOrder(self):
        return [(att, ii) for ii, att in enumerate(self._attributeNameList)]

    def getAttributeIndex(self, attributeName):
        try:
            return self._attributeNameList.index(attributeName)
        except Exception:
            return -1

    def hasAttribute(self, attributeName):
        return attributeName in self._attributeNameList

    def getIndex(self, attributeName):
        try:
            return self._attributeNameList.index(attributeName)
        except Exception:
            return -1

    def getItemNameList(self):
        return ["_" + self._name + "." + att for att in self._attributeNameList]

    def append(self, row):
        self._rowList.append(row)

    def appendAttribute(self, attributeName):
        attributeNameLC = attributeName.lower()
        if attributeNameLC in self._catalog:
            i = self._attributeNameList.index(self._catalog[attributeNameLC])
            self._attributeNameList[i] = attributeName
            self._catalog[attributeNameLC] = attributeName
        else:
            self._attributeNameList.append(attributeName)
            self._catalog[attributeNameLC] = attributeName
        self._numAttributes = len(self._attributeNameList)

    def appendAttributeExtendRows(self, attributeName):
        attributeNameLC = attributeName.lower()
        if attributeNameLC in self._catalog:
            i = self._attributeNameList.index(self._catalog[attributeNameLC])
            self._attributeNameList[i] = attributeName
            self._catalog[attributeNameLC] = attributeName
        else:
            self._attributeNameList.append(attributeName)
            self._catalog[attributeNameLC] = attributeName
            # Add placeholder to existing rows for the new attribute
            for row in self._rowList:
                row.append("?")
        self._numAttributes = len(self._attributeNameList)

    def getValue(self, attributeName=None, rowIndex=None):
        attribute = attributeName if attributeName is not None else self.__currentAttribute
        rowI = rowIndex if rowIndex is not None else self.__currentRowIndex

        if isinstance(attribute, str) and isinstance(rowI, int):
            try:
                return self._rowList[rowI][self._attributeNameList.index(attribute)]
            except IndexError:
                raise IndexError
        raise IndexError(attribute)

    def setValue(self, value, attributeName=None, rowIndex=None):
        attribute = attributeName if attributeName is not None else self.__currentAttribute
        rowI = rowIndex if rowIndex is not None else self.__currentRowIndex

        if isinstance(attribute, str) and isinstance(rowI, int):
            try:
                # Extend row list if needed
                for _ in range(rowI + 1 - len(self._rowList)):
                    self._rowList.append(self.__emptyRow())
                ll = len(self._rowList[rowI])
                ind = self._attributeNameList.index(attribute)
                # Extend the row if the attribute index is out of range
                if ind >= ll:
                    self._rowList[rowI].extend([None] * (ind - ll + 1))
                self._rowList[rowI][ind] = value
            except IndexError:
                self.__lfh.write(
                    "DataCategory(setValue) index error category %s attribute %s "
                    "index %d value %r\n" % (self._name, attribute, rowI, value))
                traceback.print_exc(file=self.__lfh)
            except ValueError:
                self.__lfh.write(
                    "DataCategory(setValue) value error category %s attribute %s "
                    "index %d value %r\n" % (self._name, attribute, rowI, value))
                traceback.print_exc(file=self.__lfh)

    def __emptyRow(self):
        return [None] * len(self._attributeNameList)

    def replaceValue(self, oldValue, newValue, attributeName):
        numReplace = 0
        if attributeName not in self._attributeNameList:
            return numReplace
        ind = self._attributeNameList.index(attributeName)
        for row in self._rowList:
            if row[ind] == oldValue:
                row[ind] = newValue
                numReplace += 1
        return numReplace

    def replaceSubstring(self, oldValue, newValue, attributeName):
        ok = False
        if attributeName not in self._attributeNameList:
            return ok
        ind = self._attributeNameList.index(attributeName)
        for row in self._rowList:
            val = row[ind]
            row[ind] = val.replace(oldValue, newValue)
            if val != row[ind]:
                ok = True
        return ok

    def invokeAttributeMethod(self, attributeName, type, method, db):
        self.__currentRowIndex = 0
        self.__currentAttribute = attributeName
        self.appendAttribute(attributeName)
        currentRowIndex = self.__currentRowIndex
        ind = self._attributeNameList.index(attributeName)

        if len(self._rowList) == 0:
            row = [None] * (len(self._attributeNameList) * 2)
            row[ind] = None
            self._rowList.append(row)

        for row in self._rowList:
            ll = len(row)
            if ind >= ll:
                row.extend([None] * (2 * ind - ll + 1))
                row[ind] = None
            exec(method.getInline())
            self.__currentRowIndex += 1
            currentRowIndex = self.__currentRowIndex

    def invokeCategoryMethod(self, type, method, db):
        self.__currentRowIndex = 0
        exec(method.getInline())

    def getAttributeLengthMaximumList(self):
        mList = [0] * len(self._attributeNameList)
        for row in self._rowList:
            for indx, val in enumerate(row):
                if val is not None:
                    mList[indx] = max(mList[indx], len(str(val)))
        return mList

    def renameAttribute(self, curAttributeName, newAttributeName):
        """Change the name of an attribute in place."""
        try:
            i = self._attributeNameList.index(curAttributeName)
            self._attributeNameList[i] = newAttributeName
            del self._catalog[curAttributeName.lower()]
            self._catalog[newAttributeName.lower()] = newAttributeName
            return True
        except Exception:
            return False

    def printIt(self, fh=sys.stdout):
        fh.write("--------------------------------------------\n")
        fh.write("  Category: %s attribute list length: %d\n" %
                 (self._name, len(self._attributeNameList)))
        for at in self._attributeNameList:
            fh.write("  Category: %s attribute: %s\n" % (self._name, at))
        fh.write("  Row value list length: %d\n" % len(self._rowList))
        for row in self._rowList[:2]:
            if len(row) == len(self._attributeNameList):
                for ii, v in enumerate(row):
                    fh.write("        %30s: %s ...\n" % (
                        self._attributeNameList[ii], str(v)[:30]))
            else:
                fh.write("+WARNING - %s data length %d attribute name length %s mismatched\n" %
                         (self._name, len(row), len(self._attributeNameList)))

    def dumpIt(self, fh=sys.stdout):
        fh.write("--------------------------------------------\n")
        fh.write("  Category: %s attribute list length: %d\n" %
                 (self._name, len(self._attributeNameList)))
        for at in self._attributeNameList:
            fh.write("  Category: %s attribute: %s\n" % (self._name, at))
        fh.write("  Value list length: %d\n" % len(self._rowList))
        for row in self._rowList:
            for ii, v in enumerate(row):
                fh.write("        %30s: %s\n" % (self._attributeNameList[ii], v))

    # ─── CIF Formatting ────────────────────────────────────────────────────────

    def __formatPdbx(self, inp):
        """Format input data following PDBx quoting rules."""
        try:
            if inp is None:
                return (["?"], 'DT_NULL_VALUE')

            if isinstance(inp, int) or self.__intRe.search(str(inp)):
                return ([str(inp)], 'DT_INTEGER')

            if isinstance(inp, float) or self.__floatRe.search(str(inp)):
                return ([str(inp)], 'DT_FLOAT')

            if inp in (".", "?"):
                return ([inp], 'DT_NULL_VALUE')

            if inp == "":
                return (["."], 'DT_NULL_VALUE')

            if not self.__wsAndQuotesRe.search(inp):
                if inp.startswith("_"):
                    return (self.__doubleQuotedList(inp), 'DT_ITEM_NAME')
                else:
                    return ([str(inp)], 'DT_UNQUOTED_STRING')
            else:
                if self.__nlRe.search(inp):
                    return (self.__semiColonQuotedList(inp), 'DT_MULTI_LINE_STRING')
                else:
                    if self.__avoidEmbeddedQuoting:
                        if not self.__dqRe.search(inp) and not self.__sqWsRe.search(inp):
                            return (self.__doubleQuotedList(inp), 'DT_DOUBLE_QUOTED_STRING')
                        elif not self.__sqRe.search(inp) and not self.__dqWsRe.search(inp):
                            return (self.__singleQuotedList(inp), 'DT_SINGLE_QUOTED_STRING')
                        else:
                            return (self.__semiColonQuotedList(inp), 'DT_MULTI_LINE_STRING')
                    else:
                        if not self.__dqRe.search(inp):
                            return (self.__doubleQuotedList(inp), 'DT_DOUBLE_QUOTED_STRING')
                        elif not self.__sqRe.search(inp):
                            return (self.__singleQuotedList(inp), 'DT_SINGLE_QUOTED_STRING')
                        else:
                            return (self.__semiColonQuotedList(inp), 'DT_MULTI_LINE_STRING')
        except Exception:
            traceback.print_exc(file=self.__lfh)

    def __dataTypePdbx(self, inp):
        """Detect the PDBx data type."""
        if inp is None:
            return 'DT_NULL_VALUE'
        if isinstance(inp, int) or self.__intRe.search(str(inp)):
            return 'DT_INTEGER'
        if isinstance(inp, float) or self.__floatRe.search(str(inp)):
            return 'DT_FLOAT'
        if inp in (".", "?"):
            return 'DT_NULL_VALUE'
        if inp == "":
            return 'DT_NULL_VALUE'
        if not self.__wsAndQuotesRe.search(inp):
            if inp.startswith("_"):
                return 'DT_ITEM_NAME'
            else:
                return 'DT_UNQUOTED_STRING'
        else:
            if self.__nlRe.search(inp):
                return 'DT_MULTI_LINE_STRING'
            else:
                if self.__avoidEmbeddedQuoting:
                    if not self.__sqRe.search(inp) and not self.__dqWsRe.search(inp):
                        return 'DT_DOUBLE_QUOTED_STRING'
                    elif not self.__dqRe.search(inp) and not self.__sqWsRe.search(inp):
                        return 'DT_SINGLE_QUOTED_STRING'
                    else:
                        return 'DT_MULTI_LINE_STRING'
                else:
                    if not self.__sqRe.search(inp):
                        return 'DT_DOUBLE_QUOTED_STRING'
                    elif not self.__dqRe.search(inp):
                        return 'DT_SINGLE_QUOTED_STRING'
                    else:
                        return 'DT_MULTI_LINE_STRING'

    def __singleQuotedList(self, inp):
        return ["'", inp, "'"]

    def __doubleQuotedList(self, inp):
        return ['"', inp, '"']

    def __semiColonQuotedList(self, inp):
        l = ["\n"]
        if inp and inp[-1] == '\n':
            l.extend([";", inp, ";", "\n"])
        else:
            l.extend([";", inp, "\n", ";", "\n"])
        return l

    def getValueFormatted(self, attributeName=None, rowIndex=None):
        attribute = attributeName if attributeName is not None else self.__currentAttribute
        rowI = rowIndex if rowIndex is not None else self.__currentRowIndex

        if isinstance(attribute, str) and isinstance(rowI, int):
            try:
                flist, ftype = self.__formatPdbx(
                    self._rowList[rowI][self._attributeNameList.index(attribute)])
                return "".join(flist)
            except IndexError:
                self.__lfh.write("attributeName %s rowI %r rowdata %r\n" % (
                    attributeName, rowI, self._rowList[rowI]))
                raise IndexError
        raise TypeError(attribute)

    def getValueFormattedByIndex(self, attributeIndex, rowIndex):
        try:
            flist, ftype = self.__formatPdbx(self._rowList[rowIndex][attributeIndex])
            return "".join(flist)
        except IndexError:
            raise IndexError

    def getAttributeValueMaxLengthList(self, steps=1):
        mList = [0] * len(self._attributeNameList)
        for row in self._rowList[::steps]:
            for indx in range(len(self._attributeNameList)):
                val = row[indx] if indx < len(row) else None
                if val is not None:
                    mList[indx] = max(mList[indx], len(str(val)))
        return mList

    def getFormatTypeList(self, steps=1):
        try:
            curDataTypeList = ['DT_NULL_VALUE'] * len(self._attributeNameList)
            for row in self._rowList[::steps]:
                for indx in range(len(self._attributeNameList)):
                    val = row[indx] if indx < len(row) else None
                    dType = self.__dataTypePdbx(val)
                    dIndx = self.__dataTypeList.index(dType)
                    cType = curDataTypeList[indx]
                    cIndx = self.__dataTypeList.index(cType)
                    cIndx = max(cIndx, dIndx)
                    curDataTypeList[indx] = self.__dataTypeList[cIndx]

            curFormatTypeList = []
            for dt in curDataTypeList:
                ii = self.__dataTypeList.index(dt)
                curFormatTypeList.append(self.__formatTypeList[ii])
        except Exception:
            self.__lfh.write(
                "PdbxDataCategory(getFormatTypeList) ++Index error at index %d in row %r\n" %
                (indx, row))
            curFormatTypeList = ['FT_UNQUOTED_STRING'] * len(self._attributeNameList)
            curDataTypeList = ['DT_UNQUOTED_STRING'] * len(self._attributeNameList)

        return curFormatTypeList, curDataTypeList

    def getFormatTypeListX(self):
        curDataTypeList = ['DT_NULL_VALUE'] * len(self._attributeNameList)
        for row in self._rowList:
            for indx in range(len(self._attributeNameList)):
                val = row[indx] if indx < len(row) else None
                dType = self.__dataTypePdbx(val)
                dIndx = self.__dataTypeList.index(dType)
                cType = curDataTypeList[indx]
                cIndx = self.__dataTypeList.index(cType)
                cIndx = max(cIndx, dIndx)
                curDataTypeList[indx] = self.__dataTypeList[cIndx]

        curFormatTypeList = []
        for dt in curDataTypeList:
            ii = self.__dataTypeList.index(dt)
            curFormatTypeList.append(self.__formatTypeList[ii])

        return curFormatTypeList, curDataTypeList
