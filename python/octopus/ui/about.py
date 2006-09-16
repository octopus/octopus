#! @PYTHON@
# -*- coding: utf-8 mode: python -*-
# about.py - About window for octopus-gui.
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

"""About window for octopus-gui.

"""

package_name=""
package_version=""

import gtk.gdk, string

import octopus.paths

class AboutWindow(gtk.AboutDialog):

    def __init__(self):
        gtk.AboutDialog.__init__(self)
        self.set_name(package_name)
        self.set_version(package_version)
        self.set_license(file("COPYING", "r").read())
        self.set_authors(map(string.strip, file("AUTHORS", "r").readlines()))
        self.set_copyright("© 2002-2006 The octopus development team")
        self.set_comments("A real-space, real-time implementation of TDDFT")
        self.set_website("http://www.tddft.org/programs/octopus/")
        self.set_website_label("Octopus Homepage")
        self.set_logo(gtk.gdk.pixbuf_new_from_file(octopus.paths.PKGDATADIR + "/oct_main.jpg"))
