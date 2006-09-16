#! @PYTHON@
# -*- coding: utf-8 mode: python -*-
# widgets.py - Widgets and other GUI elements for Octopus.
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
import gnome.ui

import about
import vartreeview

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

    def __init__(self, appname="appname", title="title"):

        # Create the toplevel window
        gnome.ui.App.__init__(self, appname, title)

        self.set_default_size(800, 600)
        self.connect('destroy', lambda w: gtk.main_quit())
        self.connect('delete_event', self.quit_delete)

        # Create a StatusBar
        import gobject
        statusbar = gobject.new(gnome.ui.AppBar, has_progress=True, has_status=True,
                                interactivity=False)
        self.statusbar=statusbar

        self.set_statusbar(statusbar)
        # Create a UIManager instance
        uimanager = gtk.UIManager()

        def select_menu(menuitem, tooltip):
            statusbar.push(tooltip)

        def deselect_menu(menuitem):
            statusbar.pop()

        def connect_menus(uimanager, action, widget):
            if isinstance(widget, gtk.MenuItem):
                tooltip = action.get_property('tooltip')
                if tooltip:
                    cid = widget.connect('select', select_menu, tooltip)
                    cid2 = widget.connect('deselect', deselect_menu)
                    widget.set_data('kiwiapp::cids', (cid, cid2))

        def disconnect_menus(uimanager, action, widget):
            if isinstance(widget, gtk.MenuItem):
                tooltip = action.get_property('tooltip')
                if tooltip:
                    cids = widget.get_data('kiwiapp::cids') or ()
                    for name, cid in cids:
                        widget.disconnect(cid)

        uimanager.connect("connect-proxy", connect_menus)
        uimanager.connect("disconnect-proxy", disconnect_menus)


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

        actiongroup0.add_toggle_actions([('ViewToolBar', None, _('View _Toolbar'), None, _('Visibility of the toolbar'), self.viewtoolbar_action, True),
                                        ('ViewStatusBar', None, _('View _Statusbar'), None, _('Visibility of the statusbar'), self.viewstatusbar_action, True)])



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

        notebook=gtk.Notebook()
        notebook.append_page(vartreeview.VariableTreeView(), gtk.Label("Variables"))
        self.set_contents(notebook)

###############

#        import gtk.gtkgl
        import vtk
        from vtk.gtk.GtkGLExtVTKRenderWindowInteractor import GtkGLExtVTKRenderWindowInteractor

        # The GtkVTKRenderWindows
        gvtk1 = GtkGLExtVTKRenderWindowInteractor()
        gvtk2 = GtkGLExtVTKRenderWindowInteractor()
        gvtk3 = GtkGLExtVTKRenderWindowInteractor()
        #gvtk1.SetDesiredUpdateRate(1000)
        gvtk1.set_size_request(800, 600)
        gvtk2.set_size_request(800, 600)
        gvtk3.set_size_request(800, 600)
        #if sys.platform != 'win32':
        #   gvtk1.set_resize_mode(gtk.RESIZE_IMMEDIATE)

        notebook.append_page(gvtk1, gtk.Label("VTK Test"))
        notebook.append_page(gvtk2, gtk.Label("Na2"))
        notebook.append_page(gvtk3, gtk.Label("benzene"))

        gvtk1.show()
        gvtk1.Initialize()
        gvtk1.Start()
        # prevents 'q' from exiting the app.
        gvtk1.AddObserver("ExitEvent", lambda o,e,x=None: x)
        gvtk2.show()
        gvtk2.Initialize()
        gvtk2.Start()
        # prevents 'q' from exiting the app.
        gvtk2.AddObserver("ExitEvent", lambda o,e,x=None: x)
        gvtk3.show()
        gvtk3.Initialize()
        gvtk3.Start()
        # prevents 'q' from exiting the app.
        gvtk3.AddObserver("ExitEvent", lambda o,e,x=None: x)

        # The VTK stuff.
        cone = vtk.vtkConeSource()
        cone.SetResolution(80)
        coneMapper = vtk.vtkPolyDataMapper()
        coneMapper.SetInput(cone.GetOutput())
        #coneActor = vtk.vtkLODActor()
        coneActor = vtk.vtkActor()
        coneActor.SetMapper(coneMapper)
        coneActor.GetProperty().SetColor(0.5, 0.5, 1.0)
        ren1 = vtk.vtkRenderer()
        gvtk1.GetRenderWindow().AddRenderer(ren1)
        ren1.AddActor(coneActor)

        import vtk, vtk.util.vtkImageImportFromArray
        import Scientific.IO.NetCDF

        densitydataNa2=Scientific.IO.NetCDF.NetCDFFile("/tmp/Na2-density-1.ncdf")
        densitydataNa2vtk=vtk.util.vtkImageImportFromArray.vtkImageImportFromArray()
        densitydataNa2vtk.SetArray(densitydataNa2.variables['rdata'].getValue())
        densityNa2Contour = vtk.vtkContourFilter()
        densityNa2Contour.SetInput(densitydataNa2vtk.GetOutput())
        densityNa2Contour.SetValue(0, 1e-7)
        densityNa2Contour.SetValue(1, 1e-4)
        densityNa2Contour.SetValue(2, 1e-3)
        densityNa2Contour.SetValue(3, 2e-3)
        densityNa2Contour.SetValue(4, 5e-3)
        densityNa2Contour.SetValue(5, 6e-3)
        densityNa2Contour.SetValue(6, 7e-3)
        densityNa2Normals = vtk.vtkPolyDataNormals()
        densityNa2Normals.SetInputConnection(densityNa2Contour.GetOutputPort())
        densityNa2Normals.SetFeatureAngle(60.0)
        densityNa2Mapper = vtk.vtkPolyDataMapper()
        densityNa2Mapper.SetInputConnection(densityNa2Normals.GetOutputPort())
        densityNa2Mapper.ScalarVisibilityOff()
        densityNa2Property = vtk.vtkProperty()
        #densityNa2Property.SetColor(colorTransferFunction)
        densityNa2Property.SetOpacity(0.30)
        densityNa2 = vtk.vtkActor()
        densityNa2.SetMapper(densityNa2Mapper)
        densityNa2.SetProperty(densityNa2Property)

        # An outline provides context around the data.
        outlineDataNa2 = vtk.vtkOutlineFilter()
        outlineDataNa2.SetInput(densitydataNa2vtk.GetOutput())
        mapOutlineNa2 = vtk.vtkPolyDataMapper()
        mapOutlineNa2.SetInputConnection(outlineDataNa2.GetOutputPort())
        outlineNa2 = vtk.vtkActor()
        outlineNa2.SetMapper(mapOutlineNa2)
        outlineNa2.GetProperty().SetColor(0, 0, 0)
        
        ren2 = vtk.vtkRenderer()
        ren2.SetBackground(1, 1, 1)
        gvtk2.GetRenderWindow().AddRenderer(ren2)
        ren2.AddActor(outlineNa2)
        ren2.AddActor(densityNa2)

        densitydatabenzene=Scientific.IO.NetCDF.NetCDFFile("/tmp/benzene-density-1.ncdf")
        densitydatabenzenevtk=vtk.util.vtkImageImportFromArray.vtkImageImportFromArray()
        densitydatabenzenevtk.SetArray(densitydatabenzene.variables['rdata'].getValue())
        densitybenzeneContour = vtk.vtkContourFilter()
        densitybenzeneContour.SetInput(densitydatabenzenevtk.GetOutput())
        densitybenzeneContour.SetValue(0, 1e-6)
        densitybenzeneContour.SetValue(1, 1e-5)
        densitybenzeneContour.SetValue(2, 1e-4)
        densitybenzeneContour.SetValue(3, 1e-3)
        densitybenzeneContour.SetValue(4, 1e-2)
        densitybenzeneContour.SetValue(5, 4e-2)
        densitybenzeneContour.SetValue(6, 0.1)
        densitybenzeneContour.SetValue(7, 0.2)
        densitybenzeneContour.SetValue(8, 0.3)
        densitybenzeneNormals = vtk.vtkPolyDataNormals()
        densitybenzeneNormals.SetInputConnection(densitybenzeneContour.GetOutputPort())
        densitybenzeneNormals.SetFeatureAngle(60.0)
        densitybenzeneMapper = vtk.vtkPolyDataMapper()
        densitybenzeneMapper.SetInputConnection(densitybenzeneNormals.GetOutputPort())
        densitybenzeneMapper.ScalarVisibilityOff()
        densitybenzeneProperty = vtk.vtkProperty()
        #densitybenzeneProperty.SetColor(colorTransferFunction)
        densitybenzeneProperty.SetOpacity(0.30)
        densitybenzene = vtk.vtkActor()
        densitybenzene.SetMapper(densitybenzeneMapper)
        densitybenzene.SetProperty(densitybenzeneProperty)

        # An outline provides context around the data.
        outlineDatabenzene = vtk.vtkOutlineFilter()
        outlineDatabenzene.SetInput(densitydatabenzenevtk.GetOutput())
        mapOutlinebenzene = vtk.vtkPolyDataMapper()
        mapOutlinebenzene.SetInputConnection(outlineDatabenzene.GetOutputPort())
        outlinebenzene = vtk.vtkActor()
        outlinebenzene.SetMapper(mapOutlinebenzene)
        outlinebenzene.GetProperty().SetColor(0, 0, 0)

        ren3 = vtk.vtkRenderer()
        ren3.SetBackground(1, 1, 1)
        gvtk3.GetRenderWindow().AddRenderer(ren3)
        ren3.AddActor(outlinebenzene)
        ren3.AddActor(densitybenzene)

###############

        self.show_all()
#        self.show()

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

    def viewtoolbar_action(self, action):
        bonobo_toolbar = self.get_dock_item_by_name('Toolbar')
        if bonobo_toolbar:
            if action.get_active():
                bonobo_toolbar.show()
            else:
                bonobo_toolbar.hide()

    def viewstatusbar_action(self, action):
        if action.get_active():
            self.statusbar.show()
        else:
            self.statusbar.hide()

    def about_action(self, action):
        aboutwindow=about.AboutWindow()
        aboutwindow.show()
        aboutwindow.run()


