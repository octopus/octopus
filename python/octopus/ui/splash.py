#! @PYTHON@
# -*- coding: utf-8 mode: python -*-
# splash.py - Splash window for octopus-gui.
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

"""Splash window for octopus-gui.

"""

import gtk
import gtk.gdk

import octopus.paths
PKGDATADIR=octopus.paths.PKGDATADIR

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
        gtk.window_set_auto_startup_notification(False)
        gtk.Window.show_now(self)
        import gobject
        gobject.timeout_add(2000, self.hide)
        while gtk.events_pending(): gtk.main_iteration()
        gtk.window_set_auto_startup_notification(True)

