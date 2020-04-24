#! /usr/bin/env python

from pylab import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
import Tkinter as tk

root=tk.Tk()
root.protocol("WM_DELETE_WINDOW", root.quit)
fig=figure(figsize=(5,4))
fig.subplots_adjust(hspace=0.30)
canvas=FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(expand=tk.TRUE, side=tk.TOP)
toolbar=NavigationToolbar2TkAgg(canvas, root)
toolbar.pack(side="bottom")

def myfunc():
   ax.cla()
   #t=arange(0, 2*pi, pi/20.)
   #ax.plot(t, sin(t))
   ax.plot(randn(100))
   canvas.draw()

bt=tk.Button(root, text="Click me", command=myfunc)
bt.pack()

ax=fig.add_subplot(111)
root.mainloop()
