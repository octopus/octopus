#! @PYTHON@
# -*- coding: utf-8 mode: python -*-
# main.py - Projects and Runs
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

"""Project and Runs

This module manages projects and runs.

"""

class OctopusRun(object):

    def __init__(self, octopusproject):
        self.name=""
        self.octopusproject=octopusproject

    

class OctopusProject(object):

    """Manage projects.

    An OctopusProject is a set of OctopusRuns that are possibly
    related to each other. But basically we only require that the
    OctopusRuns live in the same folder which is governed by this
    OctopusProject.
    
    """

    def __init__(self):
        self.runs=[]

