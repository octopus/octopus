#! @PYTHON@
# -*- coding: utf-8 mode: python -*-
# main.py - Main module for the graphical user interface
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

"""Main module for the graphical user interface.

This module manages the creation of the first gtk windows and the main
gtk loop.

"""

import gtk
import gnome.ui

import widgets

class MainApplication(object):
    """The class managing the application.

    """

    def __init__(self, package_name, package_version):
        """

        """

        gnome.init(package_name, package_version)

    def run(self):
        splashwindow=widgets.SplashScreen()
        splashwindow.show()
        projectwindow=widgets.ProjectWindow()
        projectwindow.show()

        gtk.main()

