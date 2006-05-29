#! @PYTHON@
# -*- coding: utf-8 mode: python -*-
# widgets.by - Widgets and other GUI elements for Octopus.
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

"""Widgets and other GUI elements for Octopus.

"""

import gtk
#import gtk.glade
import gtk.gdk

import gnome.ui

class SplashScreen(gtk.Window):

    def __init__(self):
        gtk.Window.__init__(self, gtk.WINDOW_TOPLEVEL)

        self.set_type_hint(gtk.gdk.WINDOW_TYPE_HINT_SPLASHSCREEN)
        self.set_position(gtk.WIN_POS_CENTER_ALWAYS)
        self.stick()
        self.set_keep_above(True)

        vbox=gtk.VBox()
#        label=gtk.Label("""octopus
#
#        I'm just an annoying spashscreen.
#        """)
#        label.show()
#        vbox.pack_start(label, True, True, 10)
        image=gtk.Image()
        image.set_from_file(PKGDATADIR + "/oct_main.jpg")
        image.show()
        vbox.pack_end(image, True, True, 10)
        self.add(vbox)
        vbox.show()

#        self.connect("button_release_event", self.hide)

    def show(self):
        gtk.Window.show_now(self)
        import gobject
        gobject.timeout_add(2000, self.hide)

class ProjectWindow(gnome.ui.App):

    ui = '''<ui>
    <menubar name="MenuBar">
      <menu action="File">
        <menuitem action="New"/>
        <menuitem action="Open"/>
        <separator />
        <menuitem action="Save"/>
        <menuitem action="SaveAs"/>
        <separator />
        <menuitem action="Print"/>
        <separator />
        <placeholder name="OldProjects"/>
        <separator />
        <menuitem action="Close"/>
        <menuitem action="Quit"/>
      </menu>
      <menu action="View">
        <menuitem action="ViewToolBar"/>
        <menuitem action="ViewStatusBar"/>
      </menu>
      <menu action="Projects">
      </menu>
      <menu action="Help">
        <menuitem action="Contents"/>
        <menuitem action="About"/>      
      </menu>
    </menubar>
    <toolbar name="Toolbar">
      <toolitem action="New"/>
      <toolitem action="Open"/>
      <toolitem action="Save"/>
      <separator />
    </toolbar>
    </ui>'''

    def __init__(self, appname="", title=""):

        # Create the toplevel window
        gnome.ui.App.__init__(self, appname, title)

        self.set_default_size(800, 600)
        self.connect('destroy', lambda w: gtk.main_quit())
        self.connect('delete_event', self.quit_delete)

        # Create a UIManager instance
        uimanager = gtk.UIManager()
        # Add the accelerator group to the toplevel window
        accelgroup = uimanager.get_accel_group()
        self.add_accel_group(accelgroup)

        # Create an ActionGroup
        actiongroup0 = gtk.ActionGroup('BaseActions')
        actiongroup1 = gtk.ActionGroup('ProjectActions')

        actiongroup0.add_actions([('File', None, _('_File')),
                                  ('Edit', None, _('_Edit')),
                                  ('View', None, _('_View')),
                                  ('Projects', None, _('_Projects')),
                                  ('Help', None, _('_Help')),
                                  ('Quit', gtk.STOCK_QUIT, _('_Quit me!'), None, _('Quit the program'), self.quit_action),
                                  ('Contents', gtk.STOCK_HELP, None, None, _('Open octopus manual'), self.noop_action),
                                  ('About', gtk.STOCK_ABOUT, None, None, _('About this application'), self.about_action)
                                  ])

        actiongroup0.add_toggle_actions([('ViewToolBar', None, _('View _ToolBar'), None, _('xxx'), self.noop_action),
                                        ('ViewStatusBar', None, _('View _StatusBar'), None, _('xxx'), self.noop_action)])



        actiongroup0.get_action('Quit').set_property('short-label', '_Quit')

        actiongroup1.add_actions([('New',  gtk.STOCK_NEW, None, None, _('Create a new project'), self.noop_action),
                                 ('Open', gtk.STOCK_OPEN, None, None, _('Open a project'), self.noop_action),
                                 ('Save', gtk.STOCK_SAVE, None, None, _('Save a project'), self.noop_action),
                                 ('SaveAs', gtk.STOCK_SAVE_AS, None, None, _('Save a project under a different name'), self.noop_action),
                                 ('Print', gtk.STOCK_PRINT, None, None, None, self.noop_action),
                                 ('Close', gtk.STOCK_CLOSE, None, None, _('Close a project'), self.noop_action),
                                 ])

        actiongroup1.get_action('Save').set_sensitive(False)
        actiongroup1.get_action('SaveAs').set_sensitive(False)
        actiongroup1.get_action('Print').set_sensitive(False)
        actiongroup1.get_action('Close').set_sensitive(False)

        # Add the actiongroup to the uimanager
        uimanager.insert_action_group(actiongroup0, -1)
        uimanager.insert_action_group(actiongroup1, -1)

        # Add a UI description
        uimanager.add_ui_from_string(self.ui)

        # Create a MenuBar
        menubar = uimanager.get_widget('/MenuBar')
        self.set_menus(menubar)

        # Create a Toolbar
        toolbar = uimanager.get_widget('/Toolbar')
        self.set_toolbar(toolbar)

        # Create a StatusBar

        import gobject
        statusbar = gobject.new(gnome.ui.AppBar, has_progress=True, has_status=True,
                                interactivity=False)
        
        self.set_statusbar(statusbar)

        self.set_contents(VariableTreeView())

#        statusbar.show()
        self.show_all()
        self.flash('sss')
#        statusbar.show()


    def noop_action(self, action):
        print _('Action not implemented yet.')
        return True

    def quit(self):
        # Checks to see if there are projects open
#        print _('Quiting? No.')
        gtk.main_quit()


    def quit_action(self, action):
        self.quit()

    def quit_delete(self, widget, event):
        self.quit()
        return True

    def about_action(self, action):
        import gtk.gdk, string
        about=gtk.AboutDialog()
        about.set_name("@PACKAGE_NAME@")
        about.set_version("@PACKAGE_VERSION@")
        about.set_license(file("COPYING", "r").read())
        about.set_authors(map(string.strip, file("AUTHORS", "r").readlines()))
        about.set_copyright("© 2002-2006 The octopus development team")
        about.set_comments("A real-space, real-time implementation of TDDFT")
        about.set_website("http://www.tddft.org/programs/octopus/")
        about.set_website_label("Octopus Homepage")
        about.set_logo(gtk.gdk.pixbuf_new_from_file(PKGDATADIR + "/oct_main.jpg"))
        about.show()
        about.run()

#    def show(self):
#        return None

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
        def MakeVars(variabletree, theroot):
            for subelement in theroot.getchildren():
                if subelement.tag=='section':
                    section=octopus.variables.base.VariableTree(subelement.get('name'))
                    MakeVars(section, subelement)
                    variabletree.add_section(section)
                elif subelement.tag=='variable':
                    variabletree.variables.append(octopus.variables.XML.element2var(subelement))
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

