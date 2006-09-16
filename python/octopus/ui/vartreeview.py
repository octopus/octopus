#! @PYTHON@
# -*- coding: utf-8 mode: python -*-
# vartreeview.py - TreeView for Octopus' variables.
# Copyright (C) 2006 A. Thimm
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
# $Id$

"""TreeView for Octopus' variables.

"""

import gtk

import octopus.paths
PKGDATADIR=octopus.paths.PKGDATADIR


def MakeVars(variabletree, theroot):
    for subelement in theroot.getchildren():
        if subelement.tag=='section':
            section=octopus.variables.base.VariableTree(subelement.get('name'))
            MakeVars(section, subelement)
            variabletree.add_section(section)
        elif subelement.tag=='variable':
            variabletree.variables.append(octopus.variables.XML.element2var(subelement))


class VariableTreeView(gtk.HPaned):

    def __init__(self):

        gtk.HPaned.__init__(self)

        # create a treestore with three columns to use as the
        # model
        self.treestore = gtk.TreeStore(str, str, object)

        import cElementTree as ET
        variablesxml=ET.parse(PKGDATADIR + "/variables.xml")
        root=variablesxml.getroot()

        import octopus.variables.base
        import octopus.variables.XML

        variabletree=octopus.variables.base.VariableTree()
        MakeVars(variabletree, root)

        def createtree2(section, treeparent):
            for variable in section.variables:
                self.treestore.append(treeparent, [variable.name, variable.type, variable])
            for subsection in section.subsections:
                createtree2(subsection, self.treestore.append(treeparent, [subsection.name, '', subsection]))
        createtree2(variabletree, None)

        # create the treeview using self.treestore
        self.treeview = gtk.TreeView(self.treestore)

        # create the Self.TreeviewColumn to display the data
        self.tvcolumn0 = gtk.TreeViewColumn('Name')
        self.tvcolumn1 = gtk.TreeViewColumn('Type')

        # add self.tvcolumn to self.treeview
        self.treeview.append_column(self.tvcolumn0)
        self.treeview.append_column(self.tvcolumn1)

        # create a CellRendererText to render the data
        cell0 = gtk.CellRendererText()
        cell1 = gtk.CellRendererText()

        # add the cell to the self.tvcolumn and allow it to expand
        self.tvcolumn0.pack_start(cell0, True)
        self.tvcolumn1.pack_start(cell1, True)

        # set the cell "text" attribute to column 0 - retrieve text
        # from that column in self.treestore
        self.tvcolumn0.add_attribute(cell0, 'text', 0)
        self.tvcolumn1.add_attribute(cell1, 'text', 1)

        # make it searchable
        self.treeview.set_search_column(0)

        # Allow sorting on the column
        self.tvcolumn0.set_sort_column_id(0)
        self.tvcolumn1.set_sort_column_id(0)

        # Disallow drag and drop reordering of rows
        self.treeview.set_reorderable(False)
#        self.treeview.set_reorderable(True)

        treeselection = self.treeview.get_selection()
        treeselection.set_mode(gtk.SELECTION_BROWSE)

        treeselection.connect("changed",self.myfunc)
        
        self.textview = gtk.TextView()
        self.textview.set_wrap_mode(gtk.WRAP_WORD)
        
        scrolled_window1 = gtk.ScrolledWindow()
        scrolled_window1.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        scrolled_window2 = gtk.ScrolledWindow()
        scrolled_window2.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        
        scrolled_window1.add_with_viewport(self.treeview)
        scrolled_window2.add_with_viewport(self.textview)
        
        self.add1(scrolled_window1)
        self.add2(scrolled_window2)

    def myfunc(self,selection):
        (model, iter) = selection.get_selected()
        row=model[iter]
        description=row[2].description
        if description!="":
            self.textview.get_buffer().set_text(row[2].description)
        else:
            self.textview.get_buffer().set_text("Missing description")
        #statusbar.push(row[0])
