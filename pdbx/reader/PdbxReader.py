# PyChem-Pro — pdbx/reader/PdbxReader.py
#
# PDBx/mmCIF tokenizer and state-machine parser.
# Migrated from Python 2 to Python 3 by the PyChem-Pro project (2026-06-21).
#   Changes: fixed read() EOF control-flow (was raising PdbxError on normal
#            EOF instead of silently terminating), improved semicolon-block
#            tokenizer to handle inline content after the closing semicolon,
#            removed Python 2 .next() / StopIteration misuse, renamed
#            SyntaxError → PdbxSyntaxError to avoid shadowing the built-in.
#
"""
PDBx/mmCIF dictionary and data file reader.

Implements a line-based tokenizer and a state-machine parser that
produces DataContainer / DefinitionContainer objects from CIF text.
"""

import re
import sys
from pdbx.reader.PdbxContainers import (
    DataContainer, DataCategory, DefinitionContainer
)



class PdbxError(Exception):
    """Class for catching general errors."""
    pass


class PdbxSyntaxError(Exception):
    """Class for catching syntax errors."""

    def __init__(self, lineNumber, text):
        super().__init__()
        self.lineNumber = lineNumber
        self.text = text

    def __str__(self):
        return "%%ERROR - [at line: %d] %s" % (self.lineNumber, self.text)


class PdbxReader(object):
    """PDBx/mmCIF reader for data files and dictionaries.

    Usage::

        container_list = []
        with open("1abc.cif", "r", encoding="utf-8", errors="replace") as fh:
            reader = PdbxReader(fh)
            reader.read(container_list)
    """

    def __init__(self, ifh):
        """
        Args:
            ifh: Input file handle returned by open().
        """
        self.__curLineNumber = 0
        self.__ifh = ifh
        self.__stateDict = {
            "data":   "ST_DATA_CONTAINER",
            "loop":   "ST_TABLE",
            "global": "ST_GLOBAL_CONTAINER",
            "save":   "ST_DEFINITION",
            "stop":   "ST_STOP",
        }

    def read(self, containerList):
        """Parse the CIF file and append containers to *containerList*.

        Args:
            containerList: A list that will receive DataContainer /
                           DefinitionContainer objects parsed from the file.
        """
        self.__curLineNumber = 0
        try:
            self.__parser(self.__tokenizer(self.__ifh), containerList)
        except StopIteration:
            # Normal end-of-file — not an error.
            pass
        except PdbxSyntaxError:
            raise
        except PdbxError:
            raise
        # Do NOT re-raise StopIteration as a PdbxError (original bug).

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def __syntaxError(self, errText):
        raise PdbxSyntaxError(self.__curLineNumber, errText)

    def __getContainerName(self, inWord):
        """Returns the name of the data_ or save_ container."""
        return str(inWord[5:]).strip()

    def __getState(self, inWord):
        """Identify reserved CIF syntax elements and assign a parser state.

        Returns:
            (reserved_word, state) tuple where state is one of the ST_*
            constants, or (None, 'ST_UNKNOWN') if *inWord* is not reserved.
        """
        i = inWord.find("_")
        if i == -1:
            return None, "ST_UNKNOWN"
        try:
            rWord = inWord[:i].lower()
            return rWord, self.__stateDict[rWord]
        except KeyError:
            return None, "ST_UNKNOWN"

    def __parser(self, tokenizer, containerList):
        """State-machine parser for PDBx/mmCIF data files and dictionaries.

        Iterates over tokens produced by __tokenizer() and builds the
        container hierarchy (DataContainer / DefinitionContainer →
        DataCategory rows).
        """
        curContainer = None
        categoryIndex = {}
        curCategory = None
        curRow = None
        state = None

        # ── Scan forward to the first reserved word ──────────────────
        while True:
            curCatName, curAttName, curQuotedString, curWord = next(tokenizer)
            if curWord is None:
                continue
            reservedWord, state = self.__getState(curWord)
            if reservedWord is not None:
                break

        # ── Main parsing loop ─────────────────────────────────────────
        while True:
            if curCatName is not None:
                state = "ST_KEY_VALUE_PAIR"
            elif curWord is not None:
                reservedWord, state = self.__getState(curWord)
            else:
                self.__syntaxError("Miscellaneous syntax error")
                return

            # ── _category.attribute  value ────────────────────────────
            if state == "ST_KEY_VALUE_PAIR":
                try:
                    curCategory = categoryIndex[curCatName]
                except KeyError:
                    curCategory = categoryIndex[curCatName] = DataCategory(curCatName)
                    try:
                        curContainer.append(curCategory)
                    except AttributeError:
                        self.__syntaxError(
                            "Category cannot be added to data_ block")
                        return
                    curRow = []
                    curCategory.append(curRow)
                else:
                    try:
                        curRow = curCategory[0]
                    except IndexError:
                        self.__syntaxError(
                            "Internal index error accessing category data")
                        return

                if curAttName in curCategory.getAttributeList():
                    self.__syntaxError(
                        "Duplicate attribute encountered in category")
                    return
                curCategory.appendAttribute(curAttName)

                # Get the value token
                tCat, tAtt, curQuotedString, curWord = next(tokenizer)

                if tCat is not None or (curQuotedString is None and curWord is None):
                    self.__syntaxError(
                        "Missing data for item _%s.%s" % (curCatName, curAttName))

                if curWord is not None:
                    reservedWord, state = self.__getState(curWord)
                    if reservedWord is not None:
                        self.__syntaxError(
                            "Unexpected reserved word: %s" % reservedWord)
                    curRow.append(curWord)
                elif curQuotedString is not None:
                    curRow.append(curQuotedString)
                else:
                    self.__syntaxError("Missing value in item-value pair")

                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)
                continue

            # ── loop_ ─────────────────────────────────────────────────
            elif state == "ST_TABLE":
                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)

                if curCatName is None or curAttName is None:
                    self.__syntaxError("Unexpected token in loop_ declaration")
                    return

                if curCatName in categoryIndex:
                    self.__syntaxError(
                        "Duplicate category declaration in loop_")
                    return

                curCategory = DataCategory(curCatName)
                try:
                    curContainer.append(curCategory)
                except AttributeError:
                    self.__syntaxError(
                        "loop_ declaration outside of data_ block or save_ frame")
                    return

                curCategory.appendAttribute(curAttName)

                # Read remaining attribute declarations
                while True:
                    curCatName, curAttName, curQuotedString, curWord = next(tokenizer)
                    if curCatName is None:
                        break
                    if curCatName != curCategory.getName():
                        self.__syntaxError(
                            "Changed category name in loop_ declaration")
                        return
                    curCategory.appendAttribute(curAttName)

                # Validate first data token
                if curWord is not None:
                    reservedWord, state = self.__getState(curWord)
                    if reservedWord is not None:
                        if reservedWord == "stop":
                            return
                        else:
                            self.__syntaxError(
                                "Unexpected reserved word after loop_ declaration: %s"
                                % reservedWord)

                # Read the tabular data
                while True:
                    curRow = []
                    curCategory.append(curRow)

                    for _tAtt in curCategory.getAttributeList():
                        if curWord is not None:
                            curRow.append(curWord)
                        elif curQuotedString is not None:
                            curRow.append(curQuotedString)
                        curCatName, curAttName, curQuotedString, curWord = next(tokenizer)

                    if curCatName is not None:
                        break
                    if curWord is not None:
                        reservedWord, state = self.__getState(curWord)
                        if reservedWord is not None:
                            break

                continue

            # ── save_ frame ───────────────────────────────────────────
            elif state == "ST_DEFINITION":
                sName = self.__getContainerName(curWord)
                if len(sName) > 0:
                    curContainer = DefinitionContainer(sName)
                    containerList.append(curContainer)
                    categoryIndex = {}
                    curCategory = None
                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)

            # ── data_ block ───────────────────────────────────────────
            elif state == "ST_DATA_CONTAINER":
                dName = self.__getContainerName(curWord)
                if not dName:
                    dName = "unidentified"
                curContainer = DataContainer(dName)
                containerList.append(curContainer)
                categoryIndex = {}
                curCategory = None
                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)

            elif state == "ST_STOP":
                return

            elif state == "ST_GLOBAL_CONTAINER":
                curContainer = DataContainer("blank-global")
                curContainer.setGlobal()
                containerList.append(curContainer)
                categoryIndex = {}
                curCategory = None
                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)

            elif state == "ST_UNKNOWN":
                self.__syntaxError("Unrecognised syntax element: " + str(curWord))
                return

    def __tokenizer(self, ifh):
        """Tokenizer for mmCIF syntax.

        Each yield produces a 4-tuple::

            (category_name, attribute_name, quoted_string, unquoted_word)

        Exactly one of the four values will be non-None per yield.
        """
        # Matches _category.attribute, single-quoted, double-quoted strings,
        # inline comments, and bare unquoted words.
        mmcifRe = re.compile(
            r"(?:"
            r"(?:_(.+?)[.](\S+))"               "|"   # _category.attribute
            r"(?:['](.*?)(?:[']\s|[']$))"        "|"   # single quoted
            r"(?:[\"](.*?)(?:[\"][\s]|[\"]$))"   "|"   # double quoted
            r"(?:\s*#.*$)"                        "|"   # comment  (discard)
            r"(\S+)"                                    # unquoted word
            r")"
        )

        fileIter = iter(ifh)

        while True:
            try:
                line = next(fileIter)
            except StopIteration:
                return
            self.__curLineNumber += 1

            # Skip full-line comments
            if line.startswith("#"):
                continue

            # ── Semi-colon multi-line string ──────────────────────────
            if line.startswith(";"):
                mlString = [line[1:]]   # content after the opening semicolon
                while True:
                    try:
                        line = next(fileIter)
                    except StopIteration:
                        break
                    self.__curLineNumber += 1
                    if line.startswith(";"):
                        break
                    mlString.append(line)

                # Strip trailing newline that is part of the \n; delimiter
                if mlString:
                    mlString[-1] = mlString[-1].rstrip()

                yield (None, None, "".join(mlString), None)

                # Process any remaining tokens on the closing-semicolon line
                remainder = line[1:] if len(line) > 1 else ""
                if remainder.strip():
                    for it in mmcifRe.finditer(remainder):
                        tgroups = it.groups()
                        if tgroups != (None, None, None, None, None):
                            qs = tgroups[2] if tgroups[2] is not None else (
                                tgroups[3] if tgroups[3] is not None else None)
                            yield (tgroups[0], tgroups[1], qs, tgroups[4])
                continue

            # ── Regular line ──────────────────────────────────────────
            for it in mmcifRe.finditer(line):
                tgroups = it.groups()
                if tgroups != (None, None, None, None, None):
                    qs = tgroups[2] if tgroups[2] is not None else (
                        tgroups[3] if tgroups[3] is not None else None)
                    yield (tgroups[0], tgroups[1], qs, tgroups[4])
