#!/usr/bin/env python3
import sys
import os
import io
from tkinter import *
from tkinter import ttk

from options import *

class EditorTab(ttk.Frame):
    def __init__(self, parent=None, file=None):
        ttk.Frame.__init__(self, parent)
        self.pack(expand=YES, fill=BOTH)

if __name__ == '__main__':
    ''' test the EditorTab class '''
    root = Tk()
    root.title(string='程序本 v0.01')

    frame=ttk.Frame(root)
    tab=EditorTab(frame)

    frame.pack(expand=YES,fill=BOTH)
    root.mainloop()

