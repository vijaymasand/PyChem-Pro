"""
Qt Compatibility Layer

Provides unified interface for both PySide6 and PyQt6.
"""

import sys

# Try to import Qt framework (support both PySide6 and PyQt6)
try:
    from PySide6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QMenuBar, QMenu, QStatusBar, QFileDialog, QMessageBox,
        QSplitter, QTabWidget, QColorDialog, QComboBox, QPushButton,
        QFrame, QLabel, QLineEdit, QCheckBox, QTreeWidget,
        QTreeWidgetItem, QTextEdit, QProgressBar, QSpinBox,
        QDoubleSpinBox, QTableWidget, QTableWidgetItem, QHeaderView,
        QScrollArea, QGridLayout, QRadioButton, QButtonGroup,
        QSizePolicy, QSlider, QGroupBox, QDockWidget, QListWidget,
        QListWidgetItem, QDialog, QInputDialog, QFormLayout,
        QGraphicsView, QGraphicsScene, QGraphicsEllipseItem, 
        QGraphicsLineItem, QGraphicsRectItem, QGraphicsPolygonItem, 
        QGraphicsPathItem, QGraphicsPixmapItem, QGraphicsSimpleTextItem, 
        QGraphicsTextItem, QGraphicsItemGroup, QGraphicsProxyWidget
    )
    from PySide6.QtCore import Qt, QThread, Signal, QObject, QSettings, QTimer, QPointF, QRect, QRectF, QCoreApplication
    from PySide6.QtGui import QAction, QActionGroup, QKeySequence, QFont, QIcon, QPixmap, QPainter, QColor, QPen, QBrush, QFontMetrics, QWheelEvent, QRadialGradient, QLinearGradient, QImage, QConicalGradient, QPainterPath, QPolygonF, QTextCursor, QPaintEvent, QPalette

    QT_FRAMEWORK = "PySide6"
    
except ImportError as e:
    print(f"PySide6 import failed: {e}")
    try:
        from PyQt6.QtWidgets import (
            QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
            QMenuBar, QMenu, QStatusBar, QFileDialog, QMessageBox,
            QSplitter, QTabWidget, QColorDialog, QComboBox, QPushButton,
            QFrame, QLabel, QLineEdit, QCheckBox, QTreeWidget,
            QTreeWidgetItem, QTextEdit, QProgressBar, QSpinBox,
            QDoubleSpinBox, QTableWidget, QTableWidgetItem, QHeaderView,
            QScrollArea, QGridLayout, QRadioButton, QButtonGroup,
            QSizePolicy, QSlider, QGroupBox, QDockWidget, QListWidget,
            QListWidgetItem, QDialog, QInputDialog, QFormLayout,
            QGraphicsView, QGraphicsScene, QGraphicsEllipseItem, 
            QGraphicsLineItem, QGraphicsRectItem, QGraphicsPolygonItem, 
            QGraphicsPathItem, QGraphicsPixmapItem, QGraphicsSimpleTextItem, 
            QGraphicsTextItem, QGraphicsItemGroup, QGraphicsProxyWidget
        )
        from PyQt6.QtCore import Qt, QThread, pyqtSignal as Signal, QObject, QSettings, QTimer, QPointF, QRect, QRectF, QCoreApplication
        from PyQt6.QtGui import QAction, QActionGroup, QKeySequence, QFont, QIcon, QPixmap, QPainter, QColor, QPen, QBrush, QFontMetrics, QWheelEvent, QRadialGradient, QLinearGradient, QImage, QConicalGradient, QPainterPath, QPolygonF, QTextCursor, QPaintEvent, QPalette

        QT_FRAMEWORK = "PyQt6"
        print("Using PyQt6 framework (PySide6 failed)")
        
    except ImportError as e2:
        print(f"PyQt6 import also failed: {e2}")
        print("Error: Neither PySide6 nor PyQt6 available")
        QT_FRAMEWORK = None
        
        # Create dummy classes if no Qt available
        class QApplication:
            def __init__(self, argv=None): pass
            def exec(self): return 0
            def quit(self): pass
            def processEvents(self): pass
            def setOverrideCursor(self, cursor): pass
            def restoreOverrideCursor(self): pass
            def setApplicationName(self, name): pass
            def setApplicationVersion(self, version): pass
            def setOrganizationName(self, name): pass
            def setFont(self, font): pass
            @staticmethod
            def setHighDpiScaleFactorRoundingPolicy(policy): pass
        class QMainWindow:
            def __init__(self): pass
            def setWindowTitle(self, title): pass
            def setGeometry(self, x, y, w, h): pass
            def setMinimumSize(self, w, h): pass
            def resize(self, w, h): pass
            def move(self, x, y): pass
            def show(self): pass
            def menuBar(self): 
                if not hasattr(self, '_menu_bar'):
                    self._menu_bar = QMenuBar()
                return self._menu_bar
            def statusBar(self): 
                if not hasattr(self, '_status_bar'):
                    self._status_bar = QStatusBar()
                return self._status_bar
            def setStatusBar(self, bar): pass
            def setCentralWidget(self, widget): pass
            def addActions(self, actions): pass
            def close(self): pass
            def raise_(self): pass
            def activateWindow(self): pass
            def setAcceptDrops(self, accept): pass
            def setStyleSheet(self, style): pass
            def setWindowIcon(self, icon): pass
            def findChild(self, type, name): return None
            def setAttribute(self, attribute, value): pass
        class QWidget:
            def __init__(self, parent=None): pass
            def setLayout(self, layout): pass
            def setWindowTitle(self, title): pass
            def setGeometry(self, x, y, w, h): pass
            def setMinimumSize(self, w, h): pass
            def resize(self, w, h): pass
            def move(self, x, y): pass
            def show(self): pass
            def menuBar(self): 
                if not hasattr(self, '_menu_bar'):
                    self._menu_bar = QMenuBar()
                return self._menu_bar
            def statusBar(self): 
                if not hasattr(self, '_status_bar'):
                    self._status_bar = QStatusBar()
                return self._status_bar
            def setCentralWidget(self, widget): pass
            def acceptDrops(self): pass
            def setAcceptDrops(self, accept): pass
            def setStyleSheet(self, style): pass
            def setWindowIcon(self, icon): pass
            def findChild(self, type, name): return None
            def setAttribute(self, attribute, value): pass
            def setMaximumHeight(self, height): pass
            def setMinimumHeight(self, height): pass
            def setFixedHeight(self, height): pass
            def setFixedWidth(self, width): pass
            def setMinimumWidth(self, width): pass
            def setContentsMargins(self, left, top, right, bottom): pass
            def setMouseTracking(self, tracking): pass
            def close(self): pass
            def raise_(self): pass
            def activateWindow(self): pass
        class QVBoxLayout:
            def __init__(self, parent=None): pass
            def addWidget(self, widget, stretch=0): pass
            def addLayout(self, layout): pass
            def addSpacing(self, spacing): pass
            def addStretch(self, stretch=0): pass
            def setContentsMargins(self, left, top, right, bottom): pass
            def setSpacing(self, spacing): pass
        class QHBoxLayout:
            def __init__(self, parent=None): pass
            def addWidget(self, widget, stretch=0): pass
            def addLayout(self, layout): pass
            def addStretch(self, stretch=0): pass
            def setContentsMargins(self, left, top, right, bottom): pass
            def setSpacing(self, spacing): pass
        class QGroupBox:
            def __init__(self, title=None, parent=None): pass
            def setTitle(self, title): pass
        
        class QSizePolicy:
            def __init__(self): pass
        
        class QSlider:
            def __init__(self, orientation=None, parent=None): 
                self._value_changed = self.Signal()
            def setRange(self, minimum, maximum): pass
            def setValue(self, value): pass
            def setStyleSheet(self, style): pass
            @property
            def valueChanged(self): return self._value_changed
            def connect(self, slot): pass
            
            class Signal:
                def connect(self, slot): pass
                def emit(self): pass
        
        class QTextCursor:
            def __init__(self): pass
        class QMenuBar:
            def __init__(self): pass
            def addMenu(self, title): 
                menu = QMenu()
                menu.title = title
                return menu
            def actions(self): return []
        class QMenu:
            def __init__(self): pass
            def addAction(self, action): pass
            def addSeparator(self): pass
            def addMenu(self, title): 
                menu = QMenu()
                menu.title = title
                return menu
            def clear(self): pass
        class QStatusBar:
            def __init__(self): pass
            def showMessage(self, message): pass
        class QFileDialog: pass
        class QMessageBox:
            def __init__(self): pass
            def setText(self, text): pass
            def setStandardButtons(self, buttons): pass
            def exec(self): return 0
            @staticmethod
            def critical(parent, title, text, buttons=0): pass
            @staticmethod
            def information(parent, title, text, buttons=0): pass
            @staticmethod
            def warning(parent, title, text, buttons=0): pass
        class QSplitter:
            def __init__(self, orientation=None): pass
            def addWidget(self, widget): pass
            def setHandleWidth(self, width): pass
            def setSizes(self, sizes): pass
            def setStretchFactor(self, index, stretch): pass
        class QTabWidget:
            def __init__(self, parent=None): pass
            def addTab(self, widget, title): pass
            def setCurrentIndex(self, index): pass
            def currentIndex(self): return 0
            def setTabPosition(self, position): pass
            def setCornerWidget(self, widget, corner): pass
            def count(self): return 0
            def removeTab(self, index): pass
            def widget(self, index): return None
            class TabPosition:
                North = 0
        class QColorDialog: pass
        class QComboBox:
            def __init__(self, parent=None): 
                self._current_index_changed = self.Signal()
            def addItems(self, items): pass
            def setCurrentText(self, text): pass
            def currentText(self): return ""
            def setStyleSheet(self, style): pass
            @property
            def currentIndexChanged(self): return self._current_index_changed
            def connect(self, slot): pass
            
            class Signal:
                def connect(self, slot): pass
                def emit(self): pass
        class QPushButton:
            def __init__(self, text=None, parent=None): 
                self._clicked = self.Signal()
                self._font = QFont()
            def setText(self, text): pass
            def setStyleSheet(self, style): pass
            def setEnabled(self, enabled): pass
            def setObjectName(self, name): pass
            def setFixedHeight(self, height): pass
            def font(self): return self._font
            def setFont(self, font): pass
            @property
            def clicked(self): return self._clicked
            def connect(self, slot): pass
            
            class Signal:
                def connect(self, slot): pass
                def emit(self): pass
        class QFrame:
            def __init__(self, parent=None): pass
            def setFrameShape(self, shape): pass
            def setFrameShadow(self, shadow): pass
            def setStyleSheet(self, style): pass
            class Shape:
                HLine = 0
                VLine = 1
            class Shadow:
                Sunken = 0
        class QLabel:
            def __init__(self, text=None, parent=None): pass
            def setText(self, text): pass
            def setWordWrap(self, wrap): pass
            def setStyleSheet(self, style): pass
            def setFont(self, font): pass
            def setObjectName(self, name): pass
            def setAlignment(self, alignment): pass
            def setFixedWidth(self, width): pass
        class QLineEdit:
            def __init__(self, parent=None): 
                self._return_pressed = self.Signal()
            def setText(self, text): pass
            def text(self): return ""
            def setFont(self, font): pass
            def setPlaceholderText(self, text): pass
            def setStyleSheet(self, style): pass
            @property
            def returnPressed(self): return self._return_pressed
            def connect(self, slot): pass
            
            class Signal:
                def connect(self, slot): pass
                def emit(self): pass
        class QCheckBox:
            def __init__(self, text=None, parent=None): 
                self._toggled = self.Signal()
            def setText(self, text): pass
            def setChecked(self, checked): pass
            def isChecked(self): return False
            def setStyleSheet(self, style): pass
            @property
            def toggled(self): return self._toggled
            def connect(self, slot): pass
            
            class Signal:
                def connect(self, slot): pass
                def emit(self): pass
        class QTreeWidget: pass
        class QTreeWidgetItem: pass
        class QTextEdit:
            def __init__(self, parent=None): pass
            def setReadOnly(self, readonly): pass
            def setMaximumHeight(self, height): pass
            def setMinimumHeight(self, height): pass
            def append(self, text): pass
            def clear(self): pass
            def toPlainText(self): return ""
            def setFont(self, font): pass
            def setStyleSheet(self, style): pass
            def setPlaceholderText(self, text): pass
        class QProgressBar:
            def __init__(self, parent=None): pass
            def setRange(self, minimum, maximum): pass
            def setValue(self, value): pass
            def setStyleSheet(self, style): pass
            def setVisible(self, visible): pass
        class QSpinBox:
            def __init__(self, parent=None): pass
            def setRange(self, minimum, maximum): pass
            def setValue(self, value): pass
            def value(self): return 0
            def setSingleStep(self, step): pass
            def setFixedWidth(self, width): pass
            def setStyleSheet(self, style): pass
        class QDoubleSpinBox: pass
        class QTableWidget:
            def __init__(self, parent=None): 
                self._horizontal_header = QHeaderView()
            def setColumnCount(self, count): pass
            def setRowCount(self, count): pass
            def setHorizontalHeaderLabels(self, labels): pass
            def setItem(self, row, column, item): pass
            def item(self, row, column): return None
            def horizontalHeader(self): return self._horizontal_header
            def setStyleSheet(self, style): pass
            def setSectionResizeMode(self, mode): pass
        class QTableWidgetItem: pass
        class QHeaderView:
            def __init__(self): pass
            def setStretchLastSection(self, stretch): pass
            def setSectionResizeMode(self, mode): pass
        class QScrollArea:
            def __init__(self): pass
            def setWidget(self, widget): pass
            def setWidgetResizable(self, resizable): pass
            def setHorizontalScrollBarPolicy(self, policy): pass
            def setStyleSheet(self, style): pass
        class QGridLayout: pass
        class QRadioButton: pass
        class QButtonGroup: pass
        
        class QDialog(QWidget):
            def __init__(self, parent=None): pass
            def exec(self): return 0
        
        class QInputDialog:
            @staticmethod
            def getDouble(parent, title, label, value=0, min=0, max=100, decimals=1, flags=0, step=1):
                return value, False
            @staticmethod
            def getText(parent, title, label, echo=0, text="", flags=0):
                return text, False
            @staticmethod
            def getInt(parent, title, label, value=0, min=0, max=100, step=1, flags=0):
                return value, False
        
        class PythonConsole:
            def __init__(self, parent=None): pass
            def set_molecule(self, molecule): pass
            def set_viewer(self, viewer_3d, viewer_2d): pass
        
        class MolViewer3D(QWidget):
            def __init__(self, parent=None): 
                super().__init__(parent)
            def set_molecule(self, molecule): pass
            def clear(self): pass
            def reset_view(self): pass
            def set_background_color(self, color): pass
            def set_sphere_scale(self, scale): pass
            def set_stick_scale(self, scale): pass
            def set_show_spheres(self, show): pass
            def set_show_sticks(self, show): pass
            def set_show_labels(self, show): pass
            def set_color_scheme(self, scheme): pass
        
        class MolViewer2D(QWidget):
            def __init__(self, parent=None): 
                super().__init__(parent)
            def set_molecule(self, molecule): pass
            def clear(self): pass
            def set_background_color(self, color): pass
            def set_bond_width(self, width): pass
            def set_atom_size(self, size): pass
            def set_color_scheme(self, scheme): pass
        
        class QGraphicsView(QWidget):
            def __init__(self, parent=None): super().__init__(parent)
            def setScene(self, scene): pass
            def setRenderHint(self, hint, on=True): pass
            def setDragMode(self, mode): pass
            def setTransformationAnchor(self, anchor): pass
            def setResizeAnchor(self, anchor): pass
            def setViewportUpdateMode(self, mode): pass
            def scale(self, sx, sy): pass
            def resetTransform(self): pass
            def centerOn(self, *args): pass
            def fitInView(self, *args, **kwargs): pass
            def itemAt(self, *args): return None
            def scene(self): return None
            class DragMode:
                NoDrag = 0
                ScrollHandDrag = 1
                RubberBandDrag = 2
            class ViewportUpdateMode:
                FullViewportUpdate = 0
                MinimalViewportUpdate = 1
                SmartViewportUpdate = 2
                NoViewportUpdate = 3
                BoundingRectViewportUpdate = 4
        
        class QGraphicsScene(QObject):
            def __init__(self, parent=None): super().__init__()
            def addItem(self, item): pass
            def removeItem(self, item): pass
            def clear(self): pass
            def setSceneRect(self, *args): pass
            def itemsBoundingRect(self): return QRectF()
            def addEllipse(self, *args): return QGraphicsEllipseItem()
            def addLine(self, *args): return QGraphicsLineItem()
            def addPath(self, *args): return QGraphicsPathItem()
            def addPixmap(self, *args): return QGraphicsPixmapItem()
            def addPolygon(self, *args): return QGraphicsPolygonItem()
            def addRect(self, *args): return QGraphicsRectItem()
            def addSimpleText(self, *args): return QGraphicsSimpleTextItem()
            def addText(self, *args): return QGraphicsTextItem()
            def addWidget(self, *args): return QGraphicsProxyWidget()
            def createItemGroup(self, items): return QGraphicsItemGroup()

        class QGraphicsItem:
            def __init__(self, parent=None): pass
            def setPos(self, *args): pass
            def setZValue(self, z): pass
            def setVisible(self, visible): pass
            def setOpacity(self, opacity): pass
            def setRotation(self, angle): pass
            def setScale(self, scale): pass
            def setToolTip(self, text): pass
            def setCursor(self, cursor): pass
            def setAcceptHoverEvents(self, accept): pass
            def setFlag(self, flag, enabled=True): pass
            def boundingRect(self): return QRectF()
            def paint(self, painter, option, widget): pass
            class GraphicsItemFlag:
                ItemIsMovable = 1
                ItemIsSelectable = 2
                ItemIsFocusable = 4
        
        class QGraphicsEllipseItem(QGraphicsItem): pass
        class QGraphicsLineItem(QGraphicsItem): pass
        class QGraphicsPathItem(QGraphicsItem): pass
        class QGraphicsPixmapItem(QGraphicsItem): pass
        class QGraphicsPolygonItem(QGraphicsItem): pass
        class QGraphicsRectItem(QGraphicsItem): pass
        class QGraphicsSimpleTextItem(QGraphicsItem): pass
        class QGraphicsTextItem(QGraphicsItem): 
            def setHtml(self, html): pass
            def setPlainText(self, text): pass
            def setDefaultTextColor(self, color): pass
            def setFont(self, font): pass
        class QGraphicsItemGroup(QGraphicsItem): pass
        class QGraphicsProxyWidget(QGraphicsItem): pass
        
        # Qt constants (moved from duplicate Qt class)
        Horizontal = 1
        Vertical = 2
        Left = 1
        Right = 2
        Top = 1
        Bottom = 2
        HCenter = 4
        VCenter = 8
        Center = HCenter | VCenter
        AlignLeft = Left
        AlignRight = Right
        AlignTop = Top
        AlignBottom = Bottom
        AlignHCenter = HCenter
        AlignVCenter = VCenter
        AlignCenter = Center
        ItemIsEnabled = 1
        ItemIsSelectable = 2
        ItemIsEditable = 4
        ItemIsDragEnabled = 8
        ItemIsDropEnabled = 16
        ItemIsUserCheckable = 32
        ItemIsTristate = 64
        
        class WidgetAttribute:
            WA_OpaquePaintEvent = 0
            WA_TranslucentBackground = 1
            WA_NoSystemBackground = 2
        
        class Orientation:
            Horizontal = 1
            Vertical = 2
            
            class Corner:
                TopRightCorner = 1
            
            class CursorShape:
                WaitCursor = None
            
            class TabPosition:
                North = 0
            
            class KeySequence:
                StandardKey = None
                Open = None
            
            class SplitterBehavior:
                KeepSize = 0
            
            class TextInteractionFlag:
                TextEditable = 0
            
            class ItemFlag:
                ItemIsEditable = 0
                ItemIsEnabled = 1
            
            class CheckState:
                Checked = 2
                Unchecked = 0
            
            class AlignmentFlag:
                AlignLeft = 1
                AlignRight = 2
                AlignCenter = 4
            
            class SortOrder:
                AscendingOrder = 0
                DescendingOrder = 1
            
            class SizePolicy:
                Fixed = 0
                Expanding = 1
            
            class WindowModality:
                NonModal = 0
                Modal = 1
            
            class WindowType:
                Window = 1
                Dialog = 2
            
            class FrameShape:
                VLine = 5
                Shadow = None
            
            class Shadow:
                Sunken = 0
            
            class Signal:
                pass
            
            class HighDpiScaleFactorRoundingPolicy:
                PassThrough = 0
                Round = 1
                RoundPreferFloor = 2
            
            class ScrollBarPolicy:
                ScrollBarAlwaysOff = 0
                ScrollBarAlwaysOn = 1
                ScrollBarAsNeeded = 2
        
        
        class QCoreApplication:
            @staticmethod
            def processEvents(): pass
        
        class QPaintEvent:
            def __init__(self, *args): pass

        class QThread:
            def start(self): pass
            def wait(self): pass
            def isRunning(self): return False
        
        class QTimer:
            def __init__(self, parent=None): 
                self._timeout = self.Signal()
            def start(self, interval): pass
            def stop(self): pass
            @property
            def timeout(self): return self._timeout
            def singleShot(self, interval, callback): pass
            def setInterval(self, interval): pass
            
            class Signal:
                def connect(self, slot): pass
                def emit(self): pass
        
        class QFont:
            def __init__(self, family=None, pointSize=-1, weight=-1, italic=False): pass
            def setFamily(self, family): pass
            def setPointSize(self, size): pass
            def setWeight(self, weight): pass
            def setItalic(self, italic): pass
            def setStyleHint(self, hint): pass
            def setBold(self, bold): pass
            class StyleHint:
                Monospace = 0
        
        class Signal:
            def __init__(self, *args): pass
            def connect(self, slot): pass
            def emit(self, *args): pass
        
        class QObject:
            def __init__(self): pass
        
        class QSettings:
            def __init__(self, *args): pass
            def value(self, key, default=None): return default
            def setValue(self, key, value): pass
        
        class QAction:
            def __init__(self, text, parent=None): 
                self.text = text
                self.parent = parent
                self._triggered = self.Signal()
            def setText(self, text): self.text = text
            def setShortcut(self, shortcut): pass
            def setCheckable(self, checkable): pass
            def setChecked(self, checked): pass
            def setEnabled(self, enabled): pass
            @property
            def triggered(self): return self._triggered
            def connect(self, slot): pass
            def setObjectName(self, name): pass
            
            class Signal:
                def connect(self, slot): pass
                def emit(self): pass
        
        class QKeySequence:
            def __init__(self, key): pass
            def toString(self): return ""
            
            class StandardKey:
                Open = "Open"
                Save = "Save"
                New = "New"
                Close = "Close"
                Quit = "Quit"
                Help = "Help"
                About = "About"
                Preferences = "Preferences"
                Undo = "Undo"
                Redo = "Redo"
                Cut = "Cut"
                Copy = "Copy"
                Paste = "Paste"
                SelectAll = "SelectAll"
        
        class QFontMetrics:
            def __init__(self): pass
        
        class QIcon:
            def __init__(self): pass
        
        class QPixmap:
            def __init__(self): pass
        
        class QPainter:
            def __init__(self): pass
            def begin(self, widget): pass
            def end(self): pass
        
        class QColor:
            def __init__(self, *args):
                if len(args) == 1:
                    # Single argument (hex string or name)
                    pass
                elif len(args) == 3:
                    # RGB values
                    self.r, self.g, self.b = args
                elif len(args) == 4:
                    # RGBA values
                    self.r, self.g, self.b, self.a = args
        
        class QPen:
            def __init__(self): pass
        
        class QBrush:
            def __init__(self): pass
        
        class QWheelEvent:
            def __init__(self): pass
        
        class QRadialGradient:
            def __init__(self): pass
        
        class QLinearGradient:
            def __init__(self): pass
        
        class QImage:
            def __init__(self): pass
        
        class QConicalGradient:
            def __init__(self): pass
        
        class QPainterPath:
            def __init__(self): pass
        
        class QPointF:
            def __init__(self): pass
        
        class QRectF:
            def __init__(self): pass
        
        class QPolygonF:
            def __init__(self): pass
        
        class Qt:
            class WidgetAttribute:
                WA_OpaquePaintEvent = 0
                WA_TranslucentBackground = 1
                WA_NoSystemBackground = 2
            
            class Orientation:
                Horizontal = 1
                Vertical = 2
            
            class Corner:
                TopRightCorner = 1
            
            class CursorShape:
                WaitCursor = None
            
            class TabPosition:
                North = 0
            
            class KeySequence:
                StandardKey = None
                Open = None
            
            class SplitterBehavior:
                KeepSize = 0
            
            class TextInteractionFlag:
                TextEditable = 0
            
            class ItemFlag:
                ItemIsEditable = 0
                ItemIsEnabled = 1
            
            class CheckState:
                Checked = 2
                Unchecked = 0
            
            class AlignmentFlag:
                AlignLeft = 1
                AlignRight = 2
                AlignCenter = 4
            
            class SortOrder:
                AscendingOrder = 0
                DescendingOrder = 1
            
            class SizePolicy:
                Fixed = 0
                Expanding = 1
            
            class WindowModality:
                NonModal = 0
                Modal = 1
            
            class WindowType:
                Window = 1
                Dialog = 2
            
            class FrameShape:
                VLine = 5
                Shadow = None
            
            class Shadow:
                Sunken = 0
            
            class Signal:
                pass
            
            class HighDpiScaleFactorRoundingPolicy:
                PassThrough = 0
                Round = 1
                RoundPreferFloor = 2
            
            class ScrollBarPolicy:
                ScrollBarAlwaysOff = 0
                ScrollBarAlwaysOn = 1
                ScrollBarAsNeeded = 2

# Export all the Qt classes for easy import
__all__ = [
    'QApplication', 'QMainWindow', 'QWidget', 'QVBoxLayout', 'QHBoxLayout',
    'QMenuBar', 'QMenu', 'QStatusBar', 'QFileDialog', 'QMessageBox',
    'QSplitter', 'QTabWidget', 'QColorDialog', 'QComboBox', 'QPushButton',
    'QFrame', 'QLabel', 'QLineEdit', 'QCheckBox', 'QTreeWidget',
    'QTreeWidgetItem', 'QTextEdit', 'QProgressBar', 'QSpinBox',
    'QDoubleSpinBox', 'QTableWidget', 'QTableWidgetItem', 'QHeaderView',
    'QScrollArea', 'QGridLayout', 'QRadioButton', 'QButtonGroup',
    'Qt', 'QThread', 'Signal', 'QObject', 'QSettings', 'QAction', 'QActionGroup',
    'QKeySequence', 'QFont', 'QIcon', 'QPixmap', 'QPalette', 'QTimer', 'QT_FRAMEWORK',
    'QPainter', 'QColor', 'QPen', 'QBrush', 'QWheelEvent',
    'QRadialGradient', 'QLinearGradient', 'QImage', 'QConicalGradient', 
    'QPainterPath', 'QPointF', 'QRectF', 'QFontMetrics', 'QSizePolicy',
    'QSlider', 'QTextCursor', 'QPolygonF', 'QGroupBox',
    'QDockWidget', 'QListWidget', 'QListWidgetItem', 'QDialog', 'QInputDialog', 'QFormLayout',
    'QCoreApplication', 'QPaintEvent',
    'QGraphicsView', 'QGraphicsScene', 'QGraphicsEllipseItem', 
    'QGraphicsLineItem', 'QGraphicsRectItem', 'QGraphicsPolygonItem', 
    'QGraphicsPathItem', 'QGraphicsPixmapItem', 'QGraphicsSimpleTextItem', 
    'QGraphicsTextItem', 'QGraphicsItemGroup', 'QGraphicsProxyWidget',
    'PythonConsole', 'MolViewer3D', 'MolViewer2D'
]

# Deferred imports of custom classes to avoid circular dependencies
try:
    from src.features.scripting_console.ui.python_console import PythonConsole
    from src.features.visualization_3d.ui.mol_viewer_3d import MolViewer3D
    from src.features.visualization_2d.ui.mol_viewer_2d import MolViewer2D
except ImportError:
    # If imports fail (e.g. during headless analysis), the dummies defined above are used
    pass
