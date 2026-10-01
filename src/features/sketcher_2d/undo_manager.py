# -*- coding: utf-8 -*-
# Adapted from ChemCanvas (GPLv3)
import copy

class UndoManager:
    MAX_UNDO_LEVELS = 50

    def __init__(self, paper):
        self.paper = paper
        self._stack = []
        self._pos = -1
        self.save_current_state("empty paper")

    def save_current_state(self, name=''):
        if self._pos < len(self._stack) - 1:
            self._stack = self._stack[:self._pos + 1]
        
        if len(self._stack) >= self.MAX_UNDO_LEVELS:
            self._stack.pop(0)
            self._pos -= 1
            
        self._stack.append(PaperState(self.paper, name))
        self._pos += 1

    def undo(self):
        if self._pos > 0:
            self._pos -= 1
            self._stack[self._pos].restore_state()
            return True
        return False

    def redo(self):
        if self._pos < len(self._stack) - 1:
            self._pos += 1
            self._stack[self._pos].restore_state()
            return True
        return False

class PaperState:
    def __init__(self, paper, name):
        self.paper = paper
        self.name = name
        self.top_levels = list(self.paper.objects)
        self.objects_data = []
        
        # Recursively collect all objects and their state
        all_objs = self._get_all_objects(self.top_levels)
        for obj in all_objs:
            state = {}
            for attr in getattr(obj, "meta__undo_properties", []):
                state[attr] = getattr(obj, attr)
            for attr in getattr(obj, "meta__undo_copy", []):
                state[attr] = copy.copy(getattr(obj, attr))
            self.objects_data.append((obj, state))

    def _get_all_objects(self, top_levels):
        objs = []
        stack = list(top_levels)
        visited = set()
        while stack:
            obj = stack.pop()
            if obj in visited: continue
            visited.add(obj)
            objs.append(obj)
            for attr in getattr(obj, "meta__undo_children_to_record", []):
                children = getattr(obj, attr)
                stack.extend(children)
        return objs

    def restore_state(self):
        # 1. Clear all current drawings from the paper
        current_objs = self._get_all_objects(self.paper.objects)
        for obj in current_objs:
            obj.clear_drawings()
        
        # 2. Restore attributes
        for obj, state in self.objects_data:
            obj.paper = self.paper
            for attr, val in state.items():
                if attr in getattr(obj, "meta__undo_copy", []):
                    setattr(obj, attr, copy.copy(val))
                else:
                    setattr(obj, attr, val)
            
            # Sync aliased attributes (e.g., vertices↔atoms, edges↔bonds)
            for alias_from, alias_to in getattr(obj, "meta__same_objects", {}).items():
                if alias_to in state:
                    setattr(obj, alias_from, getattr(obj, alias_to))
        
        # 3. Rebuild OASA internal graph structures (_neighbors)
        for obj, state in self.objects_data:
            if hasattr(obj, 'atoms') and hasattr(obj, 'bonds'):
                # It's a Molecule — rebuild OASA neighbor mappings
                for atom in obj.atoms:
                    atom._neighbors = {}
                for bond in obj.bonds:
                    if hasattr(bond, 'atoms') and len(bond.atoms) == 2:
                        a1, a2 = bond.atoms[0], bond.atoms[1]
                        a1._neighbors[bond] = a2
                        a2._neighbors[bond] = a1
                if hasattr(obj, "_flush_cache"):
                    obj._flush_cache()
        
        # 4. Set top-level objects and redraw
        self.paper.objects = list(self.top_levels)
        
        for obj, state in self.objects_data:
            obj.draw()
            
        self.paper.redraw_dirty_objects()
